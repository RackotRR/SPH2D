#include <filesystem>
#include <csv-parser/csv.hpp>

#include "Grid.h"

using namespace sphfio;
using RR::Memory::heap_darray;

static rr_uint loadLayerFromFileMM(const std::filesystem::path& path,
	ParamsPtr params,
	heap_darray<rr_float2>& r,
	heap_darray<rr_int>& itype,
	heap_darray<rr_float2>& v,
	heap_darray<rr_float>& p,
	heap_darray<rr_float>& rho) 
{
	csv::CSVReader reader{ path.string() };

	rr_uint j = 0;
	for (const auto& row : reader) {
		if (j == params->ntotal) {
			j = 0;
			break;
		}

		r(j).x = row[NAME_VARIABLE_X].get<rr_float>();
		r(j).y = row[NAME_VARIABLE_Y].get<rr_float>();
		itype(j) = row[NAME_VARIABLE_ITYPE].get<rr_int>();

		if (params->save_velocity) {
			v(j).x = row[NAME_VARIABLE_VX].get<rr_float>();
			v(j).y = row[NAME_VARIABLE_VY].get<rr_float>();
		}

		if (params->save_pressure) {
			p(j) = row[NAME_VARIABLE_P].get<rr_float>();
		}

		if (params->save_density) {
			rho(j) = row[NAME_VARIABLE_RHO].get<rr_float>();
		}

		++j;
	}

	return j;
}

TimeLayer::TimeLayer(const ExperimentLayer& experiment_layer, ParamsPtr params) :
	r{ params->ntotal },
	itype{ params->ntotal },

	v{ params->save_velocity ? params->ntotal : 0 },
	p{ params->save_velocity ? params->ntotal : 0 },
	rho{ params->save_velocity ? params->ntotal : 0 },

	ntotal{ loadLayerFromFileMM(experiment_layer.path, params,
		r,
		itype,
		v,
		p,
		rho) },

	time{ experiment_layer.get_time() }
{
}


static rr_float getParticleVx(const TimeLayer& particles, rr_uint i) { return particles.v(i).x; }
static rr_float getParticleVy(const TimeLayer& particles, rr_uint i) { return particles.v(i).y; }
static rr_float getParticleP(const TimeLayer& particles, rr_uint i) { return particles.p(i); }
static rr_float getParticleRho(const TimeLayer& particles, rr_uint i) { return particles.rho(i); }
using ParticleVarGetter = rr_float (*)(const TimeLayer& particles, rr_uint i);

static auto& getGetterByTag(const std::string& tag) {
	static const std::unordered_map<std::string, ParticleVarGetter> dict{
		{NAME_VARIABLE_VX, getParticleVx},
		{NAME_VARIABLE_VY, getParticleVy},
		{NAME_VARIABLE_P, getParticleP},
		{NAME_VARIABLE_RHO, getParticleRho}
	};

	return dict.at(tag);
}

rr_float TimeLayer::getByTag(const std::string& value, rr_uint i) const {
	auto& getter = getGetterByTag(value);
	return getter(*this, i);
}