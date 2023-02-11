#include <vector> 
#include <forward_list>
#include <stdexcept>  

#include "GridUtils.h"
#include "GridFind.h"
#include "Kernel.h"

#include <iostream>

namespace {

	// двумерная сетка для поиска соседних частиц 
	class ParticlesGrid {
		using NeighbourBlocks = std::vector<const std::vector<rr_uint>*>; // список соседних блоков
		using Particles = std::vector<rr_uint>; // список частиц в блоке
		using Grid = std::vector<std::vector<Particles>>; // сетка блоков
	private:
		Grid grid;

		const rr_uint sizeX, sizeY;
		static constexpr bool throwIfNotValidIndex = false;

		// проверка существования блока [x][y]
		bool IsIndexValid(rr_uint x, rr_uint y) const {
			return x < sizeX && y < sizeY;
		}
	public:
		// Создаём сетку блоков
		ParticlesGrid(rr_uint sizeX, rr_uint sizeY)
			: sizeX{ sizeX }, sizeY{ sizeY }
		{
			grid = std::vector<std::vector<Particles>>(sizeX);
			for (rr_uint i{}; i < sizeX; i++) {
				grid[i] = std::vector<Particles>(sizeY);
			}
		}

		// добавляет данные в нужный блок
		void AddAt(rr_uint indexX, rr_uint indexY, rr_uint indexParticle) {
			bool isIndexValid{ IsIndexValid(indexX, indexY) };

			// если нужно бросать исключения при недействительном индексе, то бросаем
			if constexpr (throwIfNotValidIndex) {
				if (!isIndexValid) {
					throw std::runtime_error{ "Index was not valid!" };
				}
			}

			if (isIndexValid) {
				grid[indexX][indexY].push_back(indexParticle);
			}
		}

		void Clear() {
			for (rr_uint row = 0; row < sizeY; row++) {
				for (rr_uint column = 0; column < sizeX; column++) {
					grid[column][row].clear();
				}
			}
		}

		// возвращает список частиц по блокам
		// результат действителен, пока существует сетка
		const NeighbourBlocks& GetNeighbours(rr_uint indexX, rr_uint indexY) const {
			// если нужно бросать исключения при недействительном индексе, то бросаем
			if constexpr (throwIfNotValidIndex) {
				if (!IsIndexValid(indexX, indexY)) {
					throw std::runtime_error{ "Index was not valid!" };
				}
			}

			static NeighbourBlocks result;
			result.clear();

			static auto AddParticlesIfIndexValid =
				[&](rr_uint x, rr_uint y) {
				if (IsIndexValid(x, y)) { // если блок существует, то добавить в результат его список частиц
					result.push_back(&grid[x][y]);
				}
			};

			// результат состоит из частиц центрального блока и 8 соседних (если они существуют)
			AddParticlesIfIndexValid(indexX - 1, indexY);
			AddParticlesIfIndexValid(indexX - 1, indexY - 1);
			AddParticlesIfIndexValid(indexX - 1, indexY + 1);
			AddParticlesIfIndexValid(indexX, indexY - 1);
			AddParticlesIfIndexValid(indexX, indexY);
			AddParticlesIfIndexValid(indexX, indexY + 1);
			AddParticlesIfIndexValid(indexX + 1, indexY);
			AddParticlesIfIndexValid(indexX + 1, indexY - 1);
			AddParticlesIfIndexValid(indexX + 1, indexY + 1);

			return result;
		}
	};


	// skale_k depends on the smoothing kernel function
	consteval rr_float GetScaleK() {
		static_assert(Params::skf > 0 && Params::skf < 4);

		rr_float scale_k;
		switch (Params::skf)
		{
		case 1:
			scale_k = 2;
			break;
		case 2:
			scale_k = 3;
			break;
		case 3:
			scale_k = 3;
			break;
		}
		return scale_k;
	}

	rr_float BlockSize() {
		return GetScaleK() * Params::hsml;
	}
	rr_uint SizeX() {
		return static_cast<rr_uint>(ceil((Params::x_maxgeom - Params::x_mingeom) / BlockSize()));
	}
	rr_uint SizeY() {
		return static_cast<rr_uint>(ceil((Params::y_maxgeom - Params::y_mingeom) / BlockSize()));
	}
}



// calculate the smoothing function for each particle and the interaction parameters used by SPH algorithm.
// Interaction pairs are determined by constucting mesh
// comparing distance with the corresponding smoothing length within nearest blocks of particles
void grid_find(
	const rr_uint ntotal, // number of particles 
	const heap_array<rr_float2, Params::maxn>& r,	// coordinates of all particles
	const heap_array<rr_int, Params::maxn>& itype, // material type: 2 - water, 0 - doesn't exist, -2 - virtual
	rr_uint& niac, // out number of interaction pairs
	heap_array<rr_uint, Params::max_interaction>& pair_i, // out, list of first partner of interaction pair
	heap_array<rr_uint, Params::max_interaction>& pair_j, // out, list of second partner of interaction pair
	heap_array<rr_float, Params::max_interaction>& w, // out, kernel for all interaction pairs 
	heap_array<rr_float2, Params::max_interaction>& dwdx) // out, derivative of kernel with respect to x, y, z
{
	static constexpr rr_float scale_k = GetScaleK();
	static const rr_float hsml = Params::hsml;
	static const rr_float blockSize = BlockSize();

	// create grid
	static ParticlesGrid grid(SizeX(), SizeY());
	grid.Clear();

	// get block index by particle index
	auto xBlock = [&](rr_uint i) {
		return static_cast<rr_uint>((-Params::x_mingeom + r(i).x) / blockSize);
	};
	auto yBlock = [&](rr_uint i) {
		return static_cast<rr_uint>((-Params::y_mingeom + r(i).y) / blockSize);
	};

	// fill in grid
	for (rr_uint i = 0; i < ntotal; i++) {
		if (itype(i) == 0) continue; // particle doesn't exist

		// block index of particle
		rr_uint indexX = xBlock(i);
		rr_uint indexY = yBlock(i);
		grid.AddAt(indexX, indexY, i);
	}

	// grid search
	niac = 0;
	for (rr_uint i = 0; i < ntotal; i++) {
		if (itype(i) == 0) continue; // particle doesn't exist
		// block index of particle
		rr_uint indexX = xBlock(i);
		rr_uint indexY = yBlock(i);

		auto& neighbourBlocks = grid.GetNeighbours(indexX, indexY);
		for (auto* block : neighbourBlocks) {
			for (rr_uint j : *block) {
				if (i >= j) continue;

				// distance between particles i and j
				rr_float2 dij = r(i) - r(j);
				rr_float dist = length(dij);

				// are neighbours?
				if (dist < scale_k * hsml) {
					if (niac < Params::max_interaction) {
						pair_i(niac) = i;
						pair_j(niac) = j;

						// kernel and derivation of kernel
						kernel(dist, dij, w(niac), dwdx(niac));
						niac++;
					}
					else {
						throw std::runtime_error{ "Too many interactions!" };
					}
				}
			}
		}
	}
}


static void make_grid(
	const rr_uint ntotal,
	const heap_array<rr_float2, Params::maxn>& r,	// coordinates of all particles
	heap_array<rr_uint, Params::maxn>& grid,
	heap_array<rr_uint, Params::max_cells>& cells_start_in_grid) // grid index of particle
{
	static heap_array<unsigned, Params::maxn> unsorted_grid;

	unsorted_grid.fill(0);
	cells_start_in_grid.fill(0);

	for (rr_uint i = 0; i < ntotal; ++i) {
		unsigned cell_idx = get_cell_idx(r(i));
		unsorted_grid(i) = cell_idx;

		cells_start_in_grid(cell_idx)++;
	}

	for (rr_uint i = 1; i < Params::max_cells; ++i) {
		cells_start_in_grid(i) += cells_start_in_grid(i - 1ull);
	}

	for (rr_uint i = ntotal; i > 0; --i) {
		rr_uint j = i - 1;

		unsigned cell_idx = unsorted_grid(j);
		grid(cells_start_in_grid(cell_idx) - 1ull) = j;
		cells_start_in_grid(cell_idx)--;
	}
}

void grid_find2(
	const rr_uint ntotal,
	const heap_array<rr_float2, Params::maxn>& r,
	heap_array<rr_uint, Params::maxn>& neighbours_count, // size of subarray of neighbours
	heap_array_md<rr_uint, Params::max_neighbours, Params::maxn>& neighbours, // neighbours indices
	heap_array_md<rr_float, Params::max_neighbours, Params::maxn>& w, // precomputed kernel
	heap_array_md<rr_float2, Params::max_neighbours, Params::maxn>& dwdr) // precomputed kernel derivative
{
	static heap_array<rr_uint, Params::maxn> grid;
	static heap_array<rr_uint, Params::max_cells> cell_starts_in_grid;
	make_grid(ntotal, r, grid, cell_starts_in_grid);

	constexpr rr_float scale_k = get_scale_k();
	const rr_float max_dist = scale_k * Params::hsml;

#pragma omp parallel for
	for (rr_iter j = 0; j < ntotal; j++) { // run through all particles
		neighbours_count(j) = 0;
		rr_uint center_cell_idx = get_cell_idx(r(j));

		rr_uint neighbour_cells[9];
		get_neighbouring_cells(center_cell_idx, neighbour_cells);
		for (rr_uint cell_i = 0; cell_i < 9; ++cell_i) { // run through neighbouring cells
			rr_uint cell_idx = neighbour_cells[cell_i];
			if (cell_idx == Params::max_cells) continue; // invalid cell

			for (rr_uint grid_i = cell_starts_in_grid(cell_idx); // run through all particles in cell
				grid_i < cell_starts_in_grid(cell_idx + 1ull);
				++grid_i)
			{
				rr_uint i = grid(grid_i); // index of particle
				// j - current particle; i - particle near
				if (i == j) continue;

				rr_float2 diff = r(i) - r(j);
				rr_float dist = length(diff);

				if (dist < max_dist) {
					rr_uint neighbour_id = neighbours_count(j)++;
					neighbours(neighbour_id, j) = i;

					kernel(dist, diff, w(neighbour_id, j), dwdr(neighbour_id, j));
				}
			} // grid_i
		} // cell_i
	} // j (particle itself)
}
