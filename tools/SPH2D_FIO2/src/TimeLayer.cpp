#include <filesystem>
#include <csv-parser/csv.hpp>

#include "Grid.h"

using namespace sphfio;
using RR::Memory::heap_darray;

static rr_uint loadLayerFromFileMM(const std::string& filename,
	rr_uint maxn,
	heap_darray<rr_float2>& r,
	heap_darray<rr_int>& itype,
	heap_darray<rr_float2>& v,
	heap_darray<rr_float>& p,
	heap_darray<rr_float>& rho) 
{
	csv::CSVReader reader(filename);

	rr_uint j = 0;
	for (const auto& row : reader) {
		if (j == params.maxn) {
			j = 0;
			break;
		}

		// todo: come with new solution
		r(j).x = row["x"].get<rr_float>();
		r(j).y = row["y"].get<rr_float>();
		itype(j) = row["itype"].get<rr_int>();

		if (reader.index_of("vx") != csv::CSV_NOT_FOUND) {
			v(j).x = row["vx"].get<rr_float>();
		}
		if (reader.index_of("vy") != csv::CSV_NOT_FOUND) {
			v(j).y = row["vy"].get<rr_float>();
		}
		if (reader.index_of("rho") != csv::CSV_NOT_FOUND) {
			rho(j) = row["rho"].get<rr_float>();
		}
		if (reader.index_of("p") != csv::CSV_NOT_FOUND) {
			p(j) = row["p"].get<rr_float>();
		}

		++j;
	}

	return j;
}

TimeLayer::TimeLayer(const std::string& filename, rr_uint maxn) :
	r{ maxn },
	itype{ maxn },
	v{ maxn },
	p{ maxn },
	rho{ maxn },
	ntotal{ loadLayerFromFileMM(filename, maxn,
		r,
		itype,
		v,
		p,
		rho) }
{
}


// todo: come with new solution
static constexpr char X_NAME[] = "x";
static constexpr char Y_NAME[] = "y";
static constexpr char ITYPE_NAME[] = "itype";
static constexpr char VX_NAME[] = "vx";
static constexpr char VY_NAME[] = "vy";
static constexpr char P_NAME[] = "p";
static constexpr char RHO_NAME[] = "rho";

static rr_float getParticleVx(const TimeLayer& particles, rr_uint i) { return particles.v(i).x; }
static rr_float getParticleVy(const TimeLayer& particles, rr_uint i) { return particles.v(i).y; }
static rr_float getParticleP(const TimeLayer& particles, rr_uint i) { return particles.p(i); }
static rr_float getParticleRho(const TimeLayer& particles, rr_uint i) { return particles.rho(i); }
using ParticleVarGetter = rr_float (*)(const TimeLayer& particles, rr_uint i);

static auto& getGetterByTag(const std::string& tag) {
	static const std::unordered_map<std::string, ParticleVarGetter> dict{
		{VX_NAME, getParticleVx},
		{VY_NAME, getParticleVy},
		{P_NAME, getParticleP},
		{RHO_NAME, getParticleRho},
	};

	return dict.at(tag);
}

rr_float TimeLayer::getByTag(const std::string& value, rr_uint i) const {
	auto& getter = getGetterByTag(value);
	return getter(*this, i);
}