#include "SPH2D_FIO.h"

using namespace sphfio;

Square::Square(const ParamsPtr& params) :
	origin_x{ params->x_mingeom },
	origin_y{ params->y_mingeom },
	size_x{ params->x_maxgeom - params->x_mingeom },
	size_y{ params->y_maxgeom - params->y_mingeom }
{
}

bool Square::contains(rr_float x, rr_float y) const {
    if (x < origin_x || y < origin_y) {
        return false;
    }

    rr_float max_x = origin_x + size_x;
    rr_float max_y = origin_y + size_y;
    if (x > max_x || y > max_y) {
        return false;
    }

    return true;
}
bool Square::contains(rr_float2 r) const {
    return contains(r.x, r.y);
}