#include <vector> 
#include <forward_list>
#include <stdexcept>  

#include "GridFind.h"
#include "Kernel.h"

#include <iostream>

namespace {

	// двумерная сетка для поиска соседних частиц 
	class ParticlesGrid {
		using NeighbourBlocks = std::vector<const std::vector<size_t>*>; // список соседних блоков
		using Particles = std::vector<size_t>; // список частиц в блоке
		using Grid = std::vector<std::vector<Particles>>; // сетка блоков
	private:
		Grid grid;

		const size_t sizeX, sizeY;
		static constexpr bool throwIfNotValidIndex = false;

		// проверка существования блока [x][y]
		bool IsIndexValid(size_t x, size_t y) const {
			return x < sizeX && y < sizeY;
		}
	public:
		// Создаём сетку блоков
		ParticlesGrid(size_t sizeX, size_t sizeY)
			: sizeX{ sizeX }, sizeY{ sizeY }
		{
			grid = std::vector<std::vector<Particles>>(sizeX);
			for (size_t i{}; i < sizeX; i++) {
				grid[i] = std::vector<Particles>(sizeY);
			}
		}

		// добавляет данные в нужный блок
		void AddAt(size_t indexX, size_t indexY, size_t indexParticle) {
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
			for (int row = 0; row < sizeY; row++) {
				for (int column = 0; column < sizeX; column++) {
					grid[column][row].clear();
				}
			}
		}

		// возвращает список частиц по блокам
		// результат действителен, пока существует сетка
		NeighbourBlocks GetNeighbours(size_t indexX, size_t indexY) const {
			// если нужно бросать исключения при недействительном индексе, то бросаем
			if constexpr (throwIfNotValidIndex) {
				if (!IsIndexValid(indexX, indexY)) {
					throw std::runtime_error{ "Index was not valid!" };
				}
			}

			NeighbourBlocks result;

			auto AddParticlesIfIndexValid{
				[&](size_t x, size_t y) {
					if (IsIndexValid(x, y)) { // если блок существует, то добавить в результат его список частиц
						result.push_back(&grid[x][y]);
					}
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
	consteval int GetScaleK() {
		static_assert(Params::skf > 0 && Params::skf < 4);

		int scale_k;
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

	double BlockSize() {
		return GetScaleK() * Params::hsml;
	}
	size_t SizeX() {
		return static_cast<size_t>(ceil((Params::x_maxgeom - Params::x_mingeom) / BlockSize()));
	}
	size_t SizeY() {
		return static_cast<size_t>(ceil((Params::y_maxgeom - Params::y_mingeom) / BlockSize()));
	}
}



// calculate the smoothing function for each particle and the interaction parameters used by SPH algorithm.
// Interaction pairs are determined by constucting mesh
// comparing distance with the corresponding smoothing length within nearest blocks of particles
void grid_find(
	const size_t ntotal, // number of particles 
	const heap_array_md<double, Params::dim, Params::maxn>& x,	// coordinates of all particles
	heap_array<int, Params::maxn>& itype, // material type: 2 - water, 0 - doesn't exist, -2 - virtual
	size_t& niac, // out number of interaction pairs
	heap_array<size_t, Params::max_interaction>& pair_i, // out, list of first partner of interaction pair
	heap_array<size_t, Params::max_interaction>& pair_j, // out, list of second partner of interaction pair
	heap_array<double, Params::max_interaction>& w, // out, kernel for all interaction pairs 
	heap_array_md<double, Params::dim, Params::max_interaction>& dwdx) // out, derivative of kernel with respect to x, y, z 
{
	static constexpr int scale_k = GetScaleK();
	static const double hsml = Params::hsml;

	double dijSqr, dij;
	heap_array<double, Params::dim> dij_dim;
	heap_array<double, Params::dim> tmp_dwdx;

	static const double blockSize = BlockSize();

	// create grid
	static ParticlesGrid grid(SizeX(), SizeY());
	grid.Clear();

	// get block index by particle index
	auto xBlock = [&](size_t i) {
		return static_cast<size_t>((-Params::x_mingeom + x(0, i)) / blockSize);
	};
	auto yBlock = [&](size_t i) {
		return static_cast<size_t>((-Params::y_mingeom + x(1, i)) / blockSize);
	};

	// fill in grid
	for (int i = 0; i < ntotal; i++) {
		if (itype(i) == 0) continue; // particle doesn't exist

		// block index of particle
		size_t indexX = xBlock(i);
		size_t indexY = yBlock(i);
		grid.AddAt(indexX, indexY, i);
	}

	// grid search
	niac = 0;
	for (int i = 0; i < ntotal; i++) {
		if (itype(i) == 0) continue; // particle doesn't exist
		// block index of particle
		size_t indexX = xBlock(i);
		size_t indexY = yBlock(i);

		auto neighbourBlocks{ grid.GetNeighbours(indexX, indexY) };
		for (auto* block : neighbourBlocks) {
			for (size_t j : *block) {
				if (i >= j) continue;

				// distance between particles i and j
				dijSqr = 0;
				for (int d = 0; d < Params::dim; d++) {
					dij_dim(d) = x(d, i) - x(d, j);
					dijSqr += sqr(dij_dim(d));
				}
				dij = sqrt(dijSqr);

				// are neighbours?
				if (dij < scale_k * hsml) {
					if (niac < Params::max_interaction) {
						pair_i(niac) = i;
						pair_j(niac) = j;

						// kernel and derivation of kernel
						kernel(dij, dij_dim, w(niac), tmp_dwdx);

						// save derivation of kernel
						for (size_t d{}; d < Params::dim; d++) {
							dwdx(d, niac) = tmp_dwdx(d);
						}
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
