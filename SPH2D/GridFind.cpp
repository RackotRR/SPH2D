#include <vector> 
#include <forward_list>
#include <stdexcept>  

#include "GridUtils.h"
#include "GridFind.h"
#include "Kernel.h"

#include <iostream>

namespace {

	// ��������� ����� ��� ������ �������� ������ 
	class ParticlesGrid {
		using NeighbourBlocks = std::vector<const std::vector<rr_uint>*>; // ������ �������� ������
		using Particles = std::vector<rr_uint>; // ������ ������ � �����
		using Grid = std::vector<std::vector<Particles>>; // ����� ������
	private:
		Grid grid;

		const rr_uint sizeX, sizeY;
		static constexpr bool throwIfNotValidIndex = false;

		// �������� ������������� ����� [x][y]
		bool IsIndexValid(rr_uint x, rr_uint y) const {
			return x < sizeX && y < sizeY;
		}
	public:
		// ������ ����� ������
		ParticlesGrid(rr_uint sizeX, rr_uint sizeY)
			: sizeX{ sizeX }, sizeY{ sizeY }
		{
			grid = std::vector<std::vector<Particles>>(sizeX);
			for (rr_uint i{}; i < sizeX; i++) {
				grid[i] = std::vector<Particles>(sizeY);
			}
		}

		// ��������� ������ � ������ ����
		void AddAt(rr_uint indexX, rr_uint indexY, rr_uint indexParticle) {
			bool isIndexValid{ IsIndexValid(indexX, indexY) };

			// ���� ����� ������� ���������� ��� ���������������� �������, �� �������
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

		// ���������� ������ ������ �� ������
		// ��������� ������������, ���� ���������� �����
		const NeighbourBlocks& GetNeighbours(rr_uint indexX, rr_uint indexY) const {
			// ���� ����� ������� ���������� ��� ���������������� �������, �� �������
			if constexpr (throwIfNotValidIndex) {
				if (!IsIndexValid(indexX, indexY)) {
					throw std::runtime_error{ "Index was not valid!" };
				}
			}

			static NeighbourBlocks result;
			result.clear();

			static auto AddParticlesIfIndexValid =
				[&](rr_uint x, rr_uint y) {
				if (IsIndexValid(x, y)) { // ���� ���� ����������, �� �������� � ��������� ��� ������ ������
					result.push_back(&grid[x][y]);
				}
			};

			// ��������� ������� �� ������ ������������ ����� � 8 �������� (���� ��� ����������)
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
						kernel(dij, w(niac), dwdx(niac));
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

void make_grid(
	const rr_uint ntotal,
	const heap_array<rr_float2, Params::maxn>& r,	// coordinates of all particles
	heap_array<rr_uint, Params::maxn>& grid,
	heap_array<rr_uint, Params::max_cells>& cells_start_in_grid) // grid index of particle
{
	static heap_array<unsigned, Params::maxn> unsorted_grid;

	unsorted_grid.fill(0);
	cells_start_in_grid.fill(0);

	unsigned cells_x = get_cell_x_from_coordinate(Params::x_maxgeom);
	unsigned cells_y = get_cell_y_from_coordinate(Params::y_maxgeom);
	unsigned cells_count = get_cell_idx(cells_x + 1, cells_y + 1);

	
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

		auto xmax = Params::x_maxgeom;
		auto xmin = Params::x_mingeom;
		auto ymin = Params::y_mingeom;
		auto ymax = Params::y_maxgeom;
		auto xy = r(j);

		unsigned cell_idx = unsorted_grid(j);
		unsigned cell_x = get_cell_x(cell_idx);
		unsigned cell_y = get_cell_y(cell_idx);
		auto cells_cell_idx = cells_start_in_grid(cell_idx);
		auto grid_i = cells_start_in_grid(cell_idx) - 1ull;
		grid(grid_i) = j;
		cells_start_in_grid(cell_idx)--;
	}
}