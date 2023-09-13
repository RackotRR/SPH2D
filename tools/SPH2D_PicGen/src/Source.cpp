#include <iostream>
#include <string>
#include <SDL.h>
#include <filesystem>

#include <Params.h>
#include <CommonIncl.h>
#include <Output.h>
#include <Input.h>

#include "SPH2DPicGenVersion.h"

struct RGB {
    Uint8 red;
    Uint8 green;
    Uint8 blue;

    bool operator==(const RGB& other) const = default;
    bool operator!=(const RGB& other) const = default;
};

RGB fluid_particle_color = { 0, 0, 255 };
RGB boundary_particle_color = { 0, 0, 0 };
RGB wave_generator_color = { 255, 0, 0 };


Uint32 get_pixel(SDL_Surface* surface, int x, int y) {
    int bpp = surface->format->BytesPerPixel;
    /* Here p is the address to the pixel we want to retrieve */
    Uint8* p = (Uint8*)surface->pixels + y * surface->pitch + x * bpp;

    switch (bpp)
    {
    case 1:
        return *p;

    case 2:
        return *(Uint16*)p;

    case 3:
        if (SDL_BYTEORDER == SDL_BIG_ENDIAN)
            return p[0] << 16 | p[1] << 8 | p[2];
        else
            return p[0] | p[1] << 8 | p[2] << 16;

    case 4:
        return *(Uint32*)p;

    default:
        return 0;       /* shouldn't happen, but avoids warnings */
    }
}
RGB get_rgb_at(int column, int row, SDL_Surface* surface) {
    Uint32 pixel = get_pixel(surface, column, row);
    RGB color;
    SDL_GetRGB(pixel, surface->format, &color.red, &color.green, &color.blue);
    return color;
}

auto count_total_particles(SDL_Surface* surface, const std::vector<RGB>& groups_color) {
    std::vector<rr_uint> groups_count(groups_color.size());

    for (rr_int row = 0; row < surface->h; ++row) {
        for (rr_int column = 0; column < surface->w; ++column) {

            RGB color = get_rgb_at(column, row, surface);

            for (rr_int i = 0; i < groups_color.size(); ++i) {
                if (groups_color[i] == color) {
                    ++groups_count[i];
                }
            }
        }
    }
    
    return groups_count;
}
auto get_particles_positions_on_surface(SDL_Surface* surface, const std::vector<RGB>& groups_color) {

    auto groups_count = count_total_particles(surface, groups_color);
    std::vector<heap_darray<rr_int2>> groups_positions(groups_color.size());
    for (rr_int i = 0; i < groups_positions.size(); ++i) {
        groups_positions[i] = heap_darray<rr_int2>(groups_count[i]);
        groups_count[i] = 0;
    }

    for (rr_int row = 0; row < surface->h; ++row) {
        for (rr_int column = 0; column < surface->w; ++column) {
            
            RGB color = get_rgb_at(column, row, surface);
            
            for (rr_int i = 0; i < groups_color.size(); ++i) {
                if (groups_color[i] == color) {
                    auto& group = groups_positions[i];
                    group[groups_count[i]] = { column, row };
                    ++groups_count[i];
                }
            }

        }
    }

    return groups_positions;
}

// expects:
// - params.delta
// - params.x_mingeom
// - params.y_maxgeom
// fills:
// - params.y_fluid_min
// - params.y_fluid_max
// - params.x_fluid_min
// - params.x_fluid_max
// returns:
// - particles filled
rr_uint fill_in_fluid_particles(
    const heap_darray<rr_int2>& positions, 
    rr_uint i, 
    heap_darray<rr_float2>& r, 
    heap_darray<rr_int>& itype) 
{
    rr_uint start_size = i;

    for (auto& [column, row] : positions) {
        itype(i) = params.TYPE_WATER;
        r(i).x = params.x_mingeom + column * params.delta;
        r(i).y = params.y_maxgeom - row * params.delta;

        params.y_fluid_min = std::min(params.y_fluid_min, r(i).y);
        params.y_fluid_max = std::max(params.y_fluid_max, r(i).y);
        params.x_fluid_min = std::min(params.x_fluid_min, r(i).x);
        params.x_fluid_max = std::max(params.x_fluid_max, r(i).x);

        ++i;
    }

    rr_uint end_size = i;
    return end_size - start_size;
}

// expects:
// - params.delta
// - params.x_mingeom
// - params.y_maxgeom
// fills:
// - params.y_boundary_min
// - params.y_boundary_max
// - params.x_boundary_min
// - params.x_boundary_max
// returns:
// - particles filled
rr_uint fill_in_boundary_particles(
    const heap_darray<rr_int2>& positions, 
    rr_uint i, 
    heap_darray<rr_float2>& r, 
    heap_darray<rr_int>& itype) 
{
    rr_uint start_size = i;

    for (auto& [column, row] : positions) {
        itype(i) = params.TYPE_BOUNDARY;
        r(i).x = params.x_mingeom + column * params.delta;
        r(i).y = params.y_maxgeom - row * params.delta;

        params.y_boundary_min = std::min(params.y_boundary_min, r(i).y);
        params.y_boundary_max = std::max(params.y_boundary_max, r(i).y);
        params.x_boundary_min = std::min(params.x_boundary_min, r(i).x);
        params.x_boundary_max = std::max(params.x_boundary_max, r(i).x);

        ++i;
    }

    rr_uint end_size = i;
    return end_size - start_size;
}

// fill in: 
// params.ntotal
// params.nfluid
// params.nvirt
// params.maxn
// params.x_maxgeom
// params.y_maxgeom
// params.x_fluid_min
// params.x_fluid_max
// params.y_fluid_min
// params.y_fluid_max
// params.x_boundary_min
// params.x_boundary_max
// params.y_boundary_min
// params.y_boundary_max
// params.nwm_particles_start
// params.nwm_particles_end
void generate_particles_data(rr_float x_mingeom, rr_float y_mingeom, rr_float particles_delta, const std::string& picture_path) {
    if (SDL_Init(SDL_INIT_EVERYTHING) != 0) {
        SDL_Log("Unable to initialize SDL: %s", SDL_GetError());
        return;
    }

    SDL_Surface* surface = SDL_LoadBMP(picture_path.c_str());
    if (!surface) {
        std::cerr << "Can't load texture file: " << picture_path << std::endl;
        SDL_Quit();
        return;
    }
    SDL_LockSurface(surface);
    
    std::vector<RGB> groups_colors = { fluid_particle_color, boundary_particle_color, wave_generator_color };
    auto groups_positions = get_particles_positions_on_surface(surface, groups_colors);
    rr_uint ntotal = 0;
    for (auto& group : groups_positions) {
        ntotal += group.size();
    }

    params.ntotal = ntotal;
    params.maxn = 1 << (1 + intlog2(params.ntotal));

    params.x_maxgeom = x_mingeom + surface->w * particles_delta;
    params.y_maxgeom = y_mingeom + surface->h * particles_delta;

    heap_darray<rr_float2> r(ntotal);
    heap_darray<rr_int> itype(ntotal);
    heap_darray<rr_float2> v(ntotal);
    heap_darray<rr_float> rho(ntotal, 1000);
    heap_darray<rr_float> p(ntotal);

    params.x_fluid_min = FLT_MAX;
    params.x_fluid_max = -FLT_MAX;
    params.y_fluid_min = FLT_MAX;
    params.y_fluid_max = -FLT_MAX;

    params.x_boundary_min = FLT_MAX;
    params.x_boundary_max = -FLT_MAX;
    params.y_boundary_min = FLT_MAX;
    params.y_boundary_max = -FLT_MAX;
    
    params.nfluid = fill_in_fluid_particles(groups_positions[0], 0, r, itype);
    params.nvirt = fill_in_boundary_particles(groups_positions[1], params.nfluid, r, itype);

    params.nwm_particles_start = params.nfluid + params.nvirt;
    params.nvirt += fill_in_boundary_particles(groups_positions[2], params.nwm_particles_start, r, itype);
    params.nwm_particles_end = params.nfluid + params.nvirt;

    try {
        setupOutput();
        dump(std::move(r),
            std::move(itype),
            std::move(v),
            std::move(rho),
            std::move(p),
            0);
    }
    catch (std::exception& ex) {
        std::cerr << "Output error: " << ex.what() << std::endl;
    }

    SDL_UnlockSurface(surface);
    SDL_FreeSurface(surface);
    SDL_Quit();
}

// expects:
// - params.dt_correction_method
// - params.CFL_coef
// - params.hsml
// - params.g
// - params.depth
// - params.eos_csqr_k
// - params.simulation_time
// - params.dt
// - params.save_step
// fills:
// - params.dt
// - params.maxtimestep
// - params.starttimestep
void fill_in_time_integration_params() {
    params.starttimestep = 0;

    if (params.dt_correction_method == 1)
    {
        params.dt = params.CFL_coef * params.hsml / (2.2f * sqrt(200 * params.g * params.depth * params.eos_csqr_k));
    }

    if (params.dt_correction_method == 2)
    {
        params.dt = 0;
        params.maxtimestep = 0;
    }
    else
    {
        params.maxtimestep = static_cast<rr_uint>(params.simulation_time / params.dt);
        if (params.maxtimestep % params.save_step != 0) { // fix last save step
            params.maxtimestep = params.maxtimestep + (params.save_step - params.maxtimestep % params.save_step);
        }
    }
}

// expects:
// - params.nwm
// - params.wave_number
// - params.wave_amp
// - params.depth
// - params.g
// fills:
// - params.freq
// - params.piston_amp
void fill_in_waves_generator_params()
{
    if (params.nwm == 2) {
        rr_float kd = params.wave_number * params.depth;
        params.freq = sqrt(params.wave_number * params.g * tanh(kd));
        params.piston_amp = params.wave_amp * 0.5f / sqr(sinh(kd)) * (sinh(kd) * cosh(kd) + kd);
    }
    else {
        params.freq = 0;
        params.piston_amp = 0;
    }
}

//params.starttimestep = 0;
// 
//params.depth;
//params.x_fluid_particles;
//params.y_fluid_particles;
// 
//params.freq
//params.piston_amp = 0;
//params.beach_x = 0;
// 
//params.max_cells = 0;
//params.maxtimestep = 0;
//params.dt
void generate_project(const std::string& project_name) {
    auto project_path = std::filesystem::path{ project_name };
    auto path_load = project_path / "SPH2DPicGenParams.json";
    auto path_write = project_path / "Params.json";
    auto picture_path = project_path / "Particles.bmp";

    if (!std::filesystem::exists(picture_path)) {
        throw std::runtime_error{ "No particles file provided: '" + picture_path.string() + "' expected" };
    }

    params.load(path_load.string());

    generate_particles_data(params.x_mingeom, params.y_mingeom, params.delta, picture_path.string());

    params.depth = params.y_fluid_max - params.y_fluid_min;
    params.x_fluid_particles = (params.x_fluid_max - params.x_fluid_min) / params.delta;
    params.y_fluid_particles = (params.y_fluid_max - params.y_fluid_min) / params.delta;

    params.max_cells = countCells(
        params.hsml,
        params.x_mingeom,
        params.y_mingeom,
        params.x_maxgeom,
        params.y_maxgeom);

    fill_in_time_integration_params();
    fill_in_waves_generator_params();

    params.makeJson(path_write.string());
}



int main(int argc, char* argv[]) {
    std::cout << "SPH2D_PicGen " << 
        SPH2D_PICGEN_VERSION_MAJOR << "." <<
        SPH2D_PICGEN_VERSION_MINOR << "." <<
        SPH2D_PICGEN_VERSION_PATCH << std::endl;

    std::string project_name;

    if (argc > 2) {
        std::cerr << "Invalid arguments num" << std::endl;
        return EXIT_FAILURE;
    }
    else if (argc == 2) {
        project_name = argv[1];
    }
    else {
        std::cout << "Enter project name: " << std::endl;
        std::cout << ">> ";
        std::getline(std::cin, project_name);
    }

    try {
        generate_project(project_name);
    }
    catch (std::exception& ex) {
        std::cerr << ex.what() << std::endl;
        std::cin.get();
    }

    return EXIT_SUCCESS;
}