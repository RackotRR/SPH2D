#include <iostream>
#include <string>
#include <SDL.h>
#include <filesystem>
#include <nlohmann/json.hpp>
#include <fstream>

#include <Params.h>
#include <ParamsIO.h>
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

static constexpr int type_boundary = -2;
static constexpr int type_non_existing = 0;
static constexpr int type_fluid = 2;

static PicGenParams pic_gen_params;
static ParticleParams particle_params;
static std::filesystem::path experiment_dir;


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

rr_uint fill_in_fluid_particles(
    const heap_darray<rr_int2>& positions, 
    rr_uint i, 
    heap_darray<rr_float2>& r, 
    heap_darray<rr_int>& itype) 
{
    rr_uint start_size = i;

    for (auto& [column, row] : positions) {
        itype(i) = type_fluid;
        r(i).x = pic_gen_params.x_mingeom + column * pic_gen_params.delta;
        r(i).y = particle_params.y_maxgeom - row * pic_gen_params.delta;

        particle_params.y_fluid_min = std::min(particle_params.y_fluid_min, r(i).y);
        particle_params.y_fluid_max = std::max(particle_params.y_fluid_max, r(i).y);
        particle_params.x_fluid_min = std::min(particle_params.x_fluid_min, r(i).x);
        particle_params.x_fluid_max = std::max(particle_params.x_fluid_max, r(i).x);

        ++i;
    }

    rr_uint end_size = i;
    return end_size - start_size;
}

rr_uint fill_in_boundary_particles(
    const heap_darray<rr_int2>& positions, 
    rr_uint i, 
    heap_darray<rr_float2>& r, 
    heap_darray<rr_int>& itype) 
{
    rr_uint start_size = i;

    for (auto& [column, row] : positions) {
        itype(i) = type_boundary;
        r(i).x = pic_gen_params.x_mingeom + column * pic_gen_params.delta;
        r(i).y = particle_params.y_maxgeom - row * pic_gen_params.delta;

        particle_params.y_boundary_min = std::min(particle_params.y_boundary_min, r(i).y);
        particle_params.y_boundary_max = std::max(particle_params.y_boundary_max, r(i).y);
        particle_params.x_boundary_min = std::min(particle_params.x_boundary_min, r(i).x);
        particle_params.x_boundary_max = std::max(particle_params.x_boundary_max, r(i).x);

        ++i;
    }

    rr_uint end_size = i;
    return end_size - start_size;
}

void generate_particles_data(
    rr_float x_mingeom, 
    rr_float y_mingeom, 
    rr_float particles_delta, 
    const std::filesystem::path& picture_path) 
{
    if (SDL_Init(SDL_INIT_EVERYTHING) != 0) {
        SDL_Log("Unable to initialize SDL: %s", SDL_GetError());
        return;
    }

    SDL_Surface* surface = SDL_LoadBMP(picture_path.string().c_str());
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

    particle_params.ntotal = ntotal;

    particle_params.x_maxgeom = x_mingeom + surface->w * particles_delta;
    particle_params.y_maxgeom = y_mingeom + surface->h * particles_delta;

    heap_darray<rr_float2> r(ntotal);
    heap_darray<rr_int> itype(ntotal);
    heap_darray<rr_float2> v(ntotal);
    heap_darray<rr_float> rho(ntotal, 1000);
    heap_darray<rr_float> p(ntotal);

    particle_params.x_fluid_min = FLT_MAX;
    particle_params.x_fluid_max = -FLT_MAX;
    particle_params.y_fluid_min = FLT_MAX;
    particle_params.y_fluid_max = -FLT_MAX;

    particle_params.x_boundary_min = FLT_MAX;
    particle_params.x_boundary_max = -FLT_MAX;
    particle_params.y_boundary_min = FLT_MAX;
    particle_params.y_boundary_max = -FLT_MAX;
    
    particle_params.nfluid = fill_in_fluid_particles(groups_positions[0], 0, r, itype);
    particle_params.nvirt = fill_in_boundary_particles(groups_positions[1], particle_params.nfluid, r, itype);

    particle_params.nwm_particles_start = particle_params.nfluid + particle_params.nvirt;
    particle_params.nvirt += fill_in_boundary_particles(groups_positions[2], particle_params.nwm_particles_start, r, itype);
    particle_params.nwm_particles_end = particle_params.nfluid + particle_params.nvirt;

    try {
        setupOutput(experiment_dir);
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

void create_default_pic_gen() {
    auto dir = std::filesystem::current_path() / "PicGenDefault";
    std::filesystem::create_directory(dir);
    auto path = dir / PicGenParams::filename;

    if (!std::filesystem::exists(path)) {
        std::ofstream stream{ path };
        nlohmann::json json;
        json["x_mingeom"] = 0.0;
        json["y_mingeom"] = 0.0;
        json["delta"] = 1.0;
        json["use_chess_order"] = false;
        json["rho0"] = 1000;
        stream << json.dump(4) << std::endl;
    }
}


static std::string getDirectoryToSearch() {
    std::cout << "Directory to search: " << std::endl;
    std::cout << "> ";
    std::string user_directory;
    std::getline(std::cin, user_directory);
    return user_directory;
}

static std::filesystem::path PicGenCLI(std::filesystem::path experiments_directory = std::filesystem::current_path()) {
    for (;;) {
        auto experiments = find_experiments(experiments_directory);
        if (experiments.empty()) {
            std::cout << "Can't find any experiment in here: " << experiments_directory << std::endl;
            experiments_directory = getDirectoryToSearch();
            std::cout << "Directory changed: " << experiments_directory << std::endl;
            continue;
        } // no experiments in directory

        auto experiment_indices = enumerate_experiments(experiments, ExperimentEnumerateCondition::pic_gen_params);

        std::cout << "Type [-1] to change directory to search." << std::endl;
        std::cout << "Type experiment number you want to load: " << std::endl;
        std::cout << "> ";
        int experiment_number;
        std::cin >> experiment_number;

        if (experiment_number == -1) {
            experiments_directory = getDirectoryToSearch();
            std::cout << "Directory changed: " << experiments_directory << std::endl;
            continue;
        }
        else if (experiment_number >= experiment_indices.size() || experiment_number < 0) {
            std::cout << "Wrong experiment number provided!" << std::endl;
            continue;
        }
        else {
            return experiments[experiment_indices[experiment_number]].dir;
        }
    }
}

void generate_experiment() {
    experiment_dir = PicGenCLI();
    auto picture_path = experiment_dir / "Particles.bmp";

    if (!std::filesystem::exists(picture_path)) {
        throw std::runtime_error{ "No particles file provided: '" + picture_path.string() + "' expected" };
    }

    pic_gen_params = load_pic_gen_params(experiment_dir);
    particle_params.rho0 = pic_gen_params.rho0;
    particle_params.delta = pic_gen_params.delta;
    particle_params.use_chess_order = pic_gen_params.use_chess_order;
    particle_params.boundary_delta = pic_gen_params.delta;
    particle_params.boundary_separation = pic_gen_params.delta;

    generate_particles_data(
        pic_gen_params.x_mingeom, 
        pic_gen_params.y_mingeom, 
        pic_gen_params.delta, 
        picture_path.string());

    particle_params.depth = particle_params.y_fluid_max - particle_params.y_fluid_min;
    particle_params.x_fluid_particles = (particle_params.x_fluid_max - particle_params.x_fluid_min) / pic_gen_params.delta;
    particle_params.y_fluid_particles = (particle_params.y_fluid_max - particle_params.y_fluid_min) / pic_gen_params.delta;

    params_make_particles_json(experiment_dir, particle_params);
}


int main(int argc, char* argv[]) {
    std::cout << "SPH2D_PicGen " << 
        SPH2D_PICGEN_VERSION_MAJOR << "." <<
        SPH2D_PICGEN_VERSION_MINOR << "." <<
        SPH2D_PICGEN_VERSION_PATCH << std::endl;

    create_default_pic_gen();

    try {
        generate_experiment();
    }
    catch (std::exception& ex) {
        std::cerr << ex.what() << std::endl;
        std::cin.get();
    }

    return EXIT_SUCCESS;
}