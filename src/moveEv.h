#ifndef __MOVEEV_H
#define __MOVEEV_H

#include "common.h"
#include "move.h"

class MoveEv {

public:
    MoveEv(movement &move_, Sint32 ev_) : move(move_), ev(ev_) { }

    movement move;
    Sint32 ev;
};

#endif
