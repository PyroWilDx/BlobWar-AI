#ifndef __MOVE_H
#define __MOVE_H
#include<SDL.h>

/** Move class
 */
struct movement {
    movement() { }
    
    movement(const movement &mv)
            : ox(mv.ox), oy(mv.oy), nx(mv.nx), ny(mv.ny) { }
    movement(const Uint8 oldx, const Uint8 oldy, const Uint8 newx, const Uint8 newy)
            : ox(oldx), oy(oldy), nx(newx), ny(newy) { }
    
    movement &operator=(const movement& mv) {
        ox=mv.ox;
        oy=mv.oy;
        nx=mv.nx;
        ny=mv.ny;
        return *this;
    }
    
    bool operator==(const movement &mv) const {
        return (ox == mv.ox && oy == mv.oy && nx == mv.nx && ny == mv.ny);
    }

    Uint8 ox;
    Uint8 oy;
    Uint8 nx;
    Uint8 ny;
};

namespace std {
    template<>
    struct hash<movement> {
        std::size_t operator()(const movement &mv) const {
            return std::hash<Uint8>()(mv.ox) ^ std::hash<Uint8>()(mv.oy) ^
                   std::hash<Uint8>()(mv.nx) ^ std::hash<Uint8>()(mv.ny);
        }
    };
}

struct point {
    point() { }

    point(const Uint8 x_, const Uint8 y_) : x(x_), y(y_) { }

    Uint8 x;
    Uint8 y;
};

#endif
