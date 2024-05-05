#include "strategy.h"
#include <chrono>
#include <fstream>
#include <algorithm>
#include <random>

// #define DEBUG_MODE

#define BASE_DEPTH 4

#define DEPTH_VAL_0 2
#define DEPTH_VAL_1 6

#define DEPTH_LIM_0 6

#define ALPHA_BETA_PAR_SEQ_LIMIT 4

void Strategy::algoNaif() {
    vector<movement> valid_moves = vector<movement>();
    computeValidMoves(_current_player, valid_moves);
    _saveBestMove(valid_moves[0]);
}

void Strategy::algoGlouton() {
    vector<movement> valid_moves = vector<movement>();
    computeValidMoves(_current_player, valid_moves);

    int bestScore = -1;

    for (movement &mv : valid_moves) {
        int currScore = 0;
        if (getMoveType(mv) == 1) currScore++;
        for (int i = std::max(0, mv.nx - 1); i < std::min(N, mv.nx + 1 + 1); i++) {
            for (int j = std::max(0, mv.ny - 1); j < std::min(N, mv.ny + 1 + 1); j++) {
                if (_blobs.get(i, j) == _other_player) {
                    currScore++;
                }
            }
        }
        if (currScore > bestScore) {
            _saveBestMove(mv);
            bestScore = currScore;
        }
    }
}

Sint32 Strategy::algoMinMaxRec(int d, Uint16 player) {
    if (d == 0) return estimateCurrentScore();
    
    vector<movement> valid_moves = vector<movement>();
    computeValidMoves(player, valid_moves);

    if (valid_moves.empty()) return estimateCurrentScore();

    vector<point> changedBlobs = vector<point>();
    if (player == _current_player) {
        Sint32 maxEval = -INFINITY_BL;

        for (movement &mv : valid_moves) {
            applyMove(_current_player, _other_player, mv, changedBlobs);
            Sint32 ev = algoMinMaxRec(d - 1, _other_player);
            revertMove(_current_player, _other_player, mv, changedBlobs);

            maxEval = std::max(maxEval, ev);
        }

        return maxEval;
    } else {
        Sint32 minEval = +INFINITY_BL;

        for (movement &mv : valid_moves) {
            applyMove(_other_player, _current_player, mv, changedBlobs);
            Sint32 ev = algoMinMaxRec(d - 1, _current_player);
            revertMove(_other_player, _current_player, mv, changedBlobs);

            minEval = std::min(minEval, ev);
        }

        return minEval;
    }
}

void Strategy::algoMinMax(int d) {
    vector<movement> valid_moves = vector<movement>();
    computeValidMoves(_current_player, valid_moves);

    vector<point> changedBlobs = vector<point>();
    Sint32 bestEval = -INFINITY_BL;

    for (movement &mv : valid_moves) {
        applyMove(_current_player, _other_player, mv, changedBlobs);
        Sint32 ev = algoMinMaxRec(d - 1, _other_player);
        revertMove(_current_player, _other_player, mv, changedBlobs);

        if (ev > bestEval) {
            _saveBestMove(mv);
            bestEval = ev;
        }
    }
}

movement Strategy::algoMinMax2(int d) {
    vector<movement> valid_moves = vector<movement>();
    computeValidMoves(_current_player, valid_moves);

    vector<point> changedBlobs = vector<point>();
    Sint32 bestEval = -INFINITY_BL;
    size_t bestI = 0;

    for (size_t i = 0; i < valid_moves.size(); i++) {
        applyMove(_current_player, _other_player, valid_moves[i], changedBlobs);
        Sint32 ev = algoMinMaxRec(d - 1, _other_player);
        revertMove(_current_player, _other_player, valid_moves[i], changedBlobs);

        if (ev > bestEval) {
            bestEval = ev;
            bestI = i;
        }
    }
    
    return valid_moves[bestI];
}

void Strategy::algoMinMax2Anytime() {
    int d = 2;
    while (d < 100) {
        movement bestMove = algoMinMax2(d);
        _saveBestMove(bestMove);
#ifdef DEBUG_MODE
        std::cout << "AlgoMinMax ================= Depth = " << d << " Done =================" << std::endl;
#endif
        d++;
    }
}

Sint32 Strategy::algoAlphaBetaRec(Sint32 alpha, Sint32 beta, int d, Uint16 player) {
    if (d == 0) return estimateCurrentScore();
    
    vector<movement> valid_moves = vector<movement>();
    computeValidMoves(player, valid_moves);

    if (valid_moves.empty()) return estimateCurrentScore();

    vector<point> changedBlobs = vector<point>();
    if (player == _current_player) {
        Sint32 maxEval = -INFINITY_BL;

        for (movement &mv : valid_moves) {
            applyMove(_current_player, _other_player, mv, changedBlobs);
            Sint32 ev = algoAlphaBetaRec(alpha, beta, d - 1, _other_player);
            revertMove(_current_player, _other_player, mv, changedBlobs);

            maxEval = std::max(maxEval, ev);

            alpha = std::max(alpha, ev);
            if (alpha >= beta) break;
        }

        return maxEval;
    } else {
        Sint32 minEval = +INFINITY_BL;

        for (movement &mv : valid_moves) {
            applyMove(_other_player, _current_player, mv, changedBlobs);
            Sint32 ev = algoAlphaBetaRec(alpha, beta, d - 1, _current_player);
            revertMove(_other_player, _current_player, mv, changedBlobs);

            minEval = std::min(minEval, ev);

            beta = std::min(beta, ev);
            if (alpha >= beta) break;
        }

        return minEval;
    }
}

void Strategy::algoAlphaBeta(int d) {
    vector<movement> valid_moves = vector<movement>();
    computeValidMoves(_current_player, valid_moves);

    vector<point> changedBlobs = vector<point>();
    Sint32 bestEval = -INFINITY_BL;
    Sint32 alpha = -INFINITY_BL;
    Sint32 beta = +INFINITY_BL;
    
    for (movement &mv : valid_moves) {
        applyMove(_current_player, _other_player, mv, changedBlobs);
        Sint32 ev = algoAlphaBetaRec(alpha, beta, d - 1, _other_player);
        revertMove(_current_player, _other_player, mv, changedBlobs);

        if (ev > bestEval) {
            _saveBestMove(mv);
            bestEval = ev;
        }
    }
}

movement Strategy::algoAlphaBeta2(int d) {
    vector<movement> valid_moves = vector<movement>();
    computeValidMoves(_current_player, valid_moves);

    vector<point> changedBlobs = vector<point>();
    Sint32 bestEval = -INFINITY_BL;
    size_t bestMoveIndex = 0;
    Sint32 alpha = -INFINITY_BL;
    Sint32 beta = +INFINITY_BL;
    
    size_t i = 0;
    for (movement &mv : valid_moves) {
        applyMove(_current_player, _other_player, mv, changedBlobs);
        Sint32 ev = algoAlphaBetaRec(alpha, beta, d - 1, _other_player);
        revertMove(_current_player, _other_player, mv, changedBlobs);

        if (ev > bestEval) {
            bestEval = ev;
            bestMoveIndex = i;
        }
        i++;
    }

    return valid_moves[bestMoveIndex];
}

void Strategy::algoAlphaBeta2Anytime() {
    int d = BASE_DEPTH;
    while (d < 100) {
        // bidiarray<Sint16> initialBlobs = bidiarray<Sint16>();
        // copyBiDiArray(_blobs, initialBlobs);
        movement bestMove = algoAlphaBeta2(d);
        _saveBestMove(bestMove);
        // _blobs = initialBlobs;
#ifdef DEBUG_MODE
        std::cout << "AlphaBeta2 ================= Depth = " << d << " Done =================" << std::endl;
#endif
        d++;
    }
}

void Strategy::applyMove(Uint16 ally, Uint16 rival, movement &mv, vector<point> &changedBlobs) {
    _blobs.set(mv.nx, mv.ny, ally);
    if (getMoveType(mv) != 1) _blobs.set(mv.ox, mv.oy, -1);

    for (int i = std::max(0, mv.nx - 1); i < std::min(N, mv.nx + 1 + 1); i++) {
        for (int j = std::max(0, mv.ny - 1); j < std::min(N, mv.ny + 1 + 1); j++) {
            if (_blobs.get(i, j) == rival) {
                _blobs.set(i, j, ally);
                changedBlobs.push_back(point(i, j));
            }
        }
    }
}

void Strategy::revertMove(Uint16 ally, Uint16 rival, movement &mv, vector<point> &changedBlobs) {
    _blobs.set(mv.nx, mv.ny, -1);
    _blobs.set(mv.ox, mv.oy, ally);

    for (point &position : changedBlobs) {
        _blobs.set(position.x, position.y, rival);
    }
    changedBlobs.clear();
}

vector<movement>& Strategy::computeValidMoves(Uint16 player, vector<movement>& valid_moves) const {
    movement mv(0, 0, 0, 0);
    for (mv.ox = 0; mv.ox < N; mv.ox++) {
        for (mv.oy = 0; mv.oy < N; mv.oy++) {
            if (_blobs.get(mv.ox, mv.oy) == player) {
                for (mv.nx = std::max(0, mv.ox - 2); mv.nx < std::min(N, mv.ox + 2 + 1); mv.nx++) {
                    for (mv.ny = std::max(0, mv.oy - 2); mv.ny < std::min(N, mv.oy + 2 + 1); mv.ny++) {
                        if (_holes.get(mv.nx, mv.ny) != 0) continue;
                        if (_blobs.get(mv.nx, mv.ny) == -1) {
                            valid_moves.push_back(movement(mv));
                        }
                    }
                }
            }
        }
    }
    return valid_moves;
}

Sint32 Strategy::estimateCurrentScore() const {
    Sint32 totalScore = 0;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (_blobs.get(i, j) == _current_player) totalScore++;
            else if (_blobs.get(i, j) == _other_player) totalScore--;
        }
    }
    return totalScore;
}

Sint32 Strategy::algoAlphaBetaSeqRec(bidiarray<Sint16> &blobs, Sint32 alpha, Sint32 beta, int d, Uint16 player) {
    if (d == 0) return estimateCurrentScore2(blobs);
    
    vector<movement> valid_moves = vector<movement>();
    computeValidMoves2(blobs, player, valid_moves);

    if (valid_moves.empty()) return estimateCurrentScore2(blobs);

    vector<point> changedBlobs = vector<point>();
    if (player == _current_player) {
        Sint32 maxEval = -INFINITY_BL;

        for (movement &mv : valid_moves) {
            applyMove2(blobs, _current_player, _other_player, mv, changedBlobs);
            Sint32 ev = algoAlphaBetaSeqRec(blobs, alpha, beta, d - 1, _other_player);
            revertMove2(blobs, _current_player, _other_player, mv, changedBlobs);

            maxEval = std::max(maxEval, ev);

            alpha = std::max(alpha, ev);
            if (alpha >= beta) break;
        }

        return maxEval;
    } else {
        Sint32 minEval = +INFINITY_BL;

        for (movement &mv : valid_moves) {
            applyMove2(blobs, _other_player, _current_player, mv, changedBlobs);
            Sint32 ev = algoAlphaBetaSeqRec(blobs, alpha, beta, d - 1, _current_player);
            revertMove2(blobs, _other_player, _current_player, mv, changedBlobs);

            minEval = std::min(minEval, ev);

            beta = std::min(beta, ev);
            if (alpha >= beta) break;
        }

        return minEval;
    }
}

Sint32 Strategy::algoAlphaBetaParRec(bidiarray<Sint16> &blobs, Sint32 alpha, Sint32 beta, int d, Uint16 player) {
    if (d == 0) return estimateCurrentScore2(blobs);

    if (d <= currLimDepth) {
        return algoAlphaBetaSeqRec(blobs, alpha, beta, d, player);
    }

    vector<movement> valid_moves = vector<movement>();
    computeValidMoves2(blobs, player, valid_moves);

    if (valid_moves.empty()) return estimateCurrentScore2(blobs);

    if (player == _current_player) {
        Sint32 maxEval = -INFINITY_BL;

        volatile bool flag = false;
        #pragma omp parallel for shared(flag)
        for (movement &mv : valid_moves) {
            if (flag) continue;

            bidiarray<Sint16> currBlobs = bidiarray<Sint16>();
            copyBiDiArray(blobs, currBlobs);
            applyMove0(currBlobs, _current_player, _other_player, mv);
            Sint32 ev = algoAlphaBetaParRec(currBlobs, alpha, beta, d - 1, _other_player);

            #pragma omp critical
            {
                maxEval = std::max(maxEval, ev);

                alpha = std::max(alpha, ev);
                if (alpha >= beta) {
                    flag = true;
                }
            }
        }

        return maxEval;
    } else {
        Sint32 minEval = +INFINITY_BL;

        volatile bool flag = false;
        #pragma omp parallel for shared(flag)
        for (movement &mv : valid_moves) {
            if (flag) continue;

            bidiarray<Sint16> currBlobs = bidiarray<Sint16>();
            copyBiDiArray(blobs, currBlobs);
            applyMove0(currBlobs, _other_player, _current_player, mv);
            Sint32 ev = algoAlphaBetaParRec(currBlobs, alpha, beta, d - 1, _current_player);

            #pragma omp critical
            {
                minEval = std::min(minEval, ev);

                beta = std::min(beta, ev);
                if (alpha >= beta) {
                    flag = true;
                }
            }
        }

        return minEval;
    }
}

void Strategy::algoAlphaBetaPar(bidiarray<Sint16> &blobs, int d) {
    vector<movement> valid_moves = vector<movement>();
    computeValidMoves2(blobs, _current_player, valid_moves);

    Sint32 bestEval = -INFINITY_BL;
    size_t minI = 0;
    Sint32 alpha = -INFINITY_BL;
    Sint32 beta = +INFINITY_BL;
    
    #pragma omp parallel for
    for (size_t i = 0; i < valid_moves.size(); i++) {
        bidiarray<Sint16> currBlobs = bidiarray<Sint16>();
        copyBiDiArray(blobs, currBlobs);
        applyMove0(currBlobs, _current_player, _other_player, valid_moves[i]);
        Sint32 ev = algoAlphaBetaParRec(currBlobs, alpha, beta, d - 1, _other_player);

        #pragma omp critical
        {
            if (ev > bestEval || (ev == bestEval && i < minI)) {
                _saveBestMove(valid_moves[i]);
                bestEval = ev;
                minI = i;
            }
        }
    }
}

void Strategy::applyMove0(bidiarray<Sint16> &blobs, Uint16 ally, Uint16 rival, movement &mv) {
    blobs.set(mv.nx, mv.ny, ally);
    if (getMoveType(mv) != 1) blobs.set(mv.ox, mv.oy, -1);

    for (int i = std::max(0, mv.nx - 1); i < std::min(N, mv.nx + 1 + 1); i++) {
        for (int j = std::max(0, mv.ny - 1); j < std::min(N, mv.ny + 1 + 1); j++) {
            if (blobs.get(i, j) == rival) {
                blobs.set(i, j, ally);
            }
        }
    }
}

void Strategy::applyMove2(bidiarray<Sint16> &blobs, Uint16 ally, Uint16 rival, movement &mv, vector<point> &changedBlobs) {
    blobs.set(mv.nx, mv.ny, ally);
    if (getMoveType(mv) != 1) blobs.set(mv.ox, mv.oy, -1);

    for (int i = std::max(0, mv.nx - 1); i < std::min(N, mv.nx + 1 + 1); i++) {
        for (int j = std::max(0, mv.ny - 1); j < std::min(N, mv.ny + 1 + 1); j++) {
            if (blobs.get(i, j) == rival) {
                blobs.set(i, j, ally);
                changedBlobs.push_back(point(i, j));
            }
        }
    }
}

void Strategy::revertMove2(bidiarray<Sint16> &blobs, Uint16 ally, Uint16 rival, movement &mv, vector<point> &changedBlobs, bool doClear) {
    blobs.set(mv.nx, mv.ny, -1);
    blobs.set(mv.ox, mv.oy, ally);

    for (point &position : changedBlobs) {
        blobs.set(position.x, position.y, rival);
    }
    if (doClear) changedBlobs.clear();
}

vector<movement> &Strategy::computeValidMoves2(bidiarray<Sint16> &blobs, Uint16 player, vector<movement> &valid_moves) const {
    movement mv(0, 0, 0, 0);
    for (mv.ox = 0; mv.ox < N; mv.ox++) {
        for (mv.oy = 0; mv.oy < N; mv.oy++) {
            if (blobs.get(mv.ox, mv.oy) == player) {
                for (mv.nx = std::max(0, mv.ox - 2); mv.nx < std::min(N, mv.ox + 2 + 1); mv.nx++) {
                    for (mv.ny = std::max(0, mv.oy - 2); mv.ny < std::min(N, mv.oy + 2 + 1); mv.ny++) {
                        if (_holes.get(mv.nx, mv.ny) != 0) continue;
                        if (blobs.get(mv.nx, mv.ny) == -1) {
                            valid_moves.push_back(movement(mv));
                        }
                    }
                }
            }
        }
    }
    return valid_moves;
}

Sint32 Strategy::estimateCurrentScore2(bidiarray<Sint16> &blobs) const {
    Sint32 totalScore = 0;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (blobs.get(i, j) == _current_player) totalScore++;
            else if (blobs.get(i, j) == _other_player) totalScore--;
        }
    }
    return totalScore;
}

movement Strategy::algoAlphaBetaPar2(bidiarray<Sint16> &blobs, int d) {
    vector<movement> valid_moves = vector<movement>();
    computeValidMoves2(blobs, _current_player, valid_moves);

    Sint32 bestEval = -INFINITY_BL;
    size_t minI = 0;
    Sint32 alpha = -INFINITY_BL;
    Sint32 beta = +INFINITY_BL;
    
    #pragma omp parallel for
    for (size_t i = 0; i < valid_moves.size(); i++) {
        bidiarray<Sint16> currBlobs = bidiarray<Sint16>();
        copyBiDiArray(blobs, currBlobs);
        applyMove0(currBlobs, _current_player, _other_player, valid_moves[i]);
        Sint32 ev = algoAlphaBetaParRec(currBlobs, alpha, beta, d - 1, _other_player);

        #pragma omp critical
        {
            if (ev > bestEval || (ev == bestEval && i < minI)) {
            // if (ev > bestEval) {
                bestEval = ev;
                minI = i;
            }
        }
    }

    return valid_moves[minI];
}

void Strategy::algoAlphaBetaPar2Anytime() {
    int d = BASE_DEPTH;
    while (d < 10) {
        if (d <= DEPTH_LIM_0) {
            currLimDepth = d - DEPTH_VAL_0;
        } else {
            currLimDepth = d - DEPTH_VAL_1;
        }
        movement bestMove = algoAlphaBetaPar2(_blobs, d);
        _saveBestMove(bestMove);
#ifdef DEBUG_MODE
        std::cout << "AlphaBetaPar2 ================= Depth = " << d << " Done =================" << std::endl;
#endif
        d++;
    }
}

void Strategy::computeValidMovesEv(bidiarray<Sint16> &blobs, Uint16 ally, Uint16 rival, vector<MoveEv> &valid_moves) {
    movement mv(0, 0, 0, 0);
    vector<point> changedBlobs = vector<point>();
    changedBlobs.reserve(8);
    for (mv.ox = 0; mv.ox < N; mv.ox++) {
        for (mv.oy = 0; mv.oy < N; mv.oy++) {
            if (blobs.get(mv.ox, mv.oy) == ally) {
                for (mv.nx = std::max(0, mv.ox - 2); mv.nx < std::min(N, mv.ox + 2 + 1); mv.nx++) {
                    for (mv.ny = std::max(0, mv.oy - 2); mv.ny < std::min(N, mv.oy + 2 + 1); mv.ny++) {
                        if (_holes.get(mv.nx, mv.ny) != 0) continue;
                        if (blobs.get(mv.nx, mv.ny) == -1) {
                            applyMove2(blobs, ally, rival, mv, changedBlobs);
                            Sint32 ev = estimateCurrentScore2(blobs);
                            valid_moves.push_back(MoveEv(mv, ev));
                            revertMove2(blobs, ally, rival, mv, changedBlobs);
                        }
                    }
                }
            }
        }
    }
}

Sint32 Strategy::algoAlphaBetaSortParRec(bidiarray<Sint16> &blobs, Sint32 alpha, Sint32 beta, int d, Uint16 player) {
    if (d == 0) return estimateCurrentScore2(blobs);

    if (d <= currLimDepth) {
        return algoAlphaBetaSeqRec(blobs, alpha, beta, d, player);
    }

    vector<MoveEv> valid_moves = vector<MoveEv>();
    computeValidMovesEv(blobs, player, (player == 1) ? 0 : 1, valid_moves);
    if (player == _current_player) {
        std::sort(valid_moves.begin(), valid_moves.end(), 
                [](const MoveEv &m1, const MoveEv &m2) { return m1.ev > m2.ev; });
    } else {
        std::sort(valid_moves.begin(), valid_moves.end(), 
                [](const MoveEv &m1, const MoveEv &m2) { return m1.ev < m2.ev; });
    }

    if (valid_moves.empty()) return estimateCurrentScore2(blobs);

    if (player == _current_player) {
        Sint32 maxEval = -INFINITY_BL;

        volatile bool flag = false;
        #pragma omp parallel for shared(flag)
        for (MoveEv &mv : valid_moves) {
            if (flag) continue;

            bidiarray<Sint16> currBlobs = bidiarray<Sint16>();
            copyBiDiArray(blobs, currBlobs);
            applyMove0(currBlobs, _current_player, _other_player, mv.move);
            Sint32 ev = algoAlphaBetaSortParRec(currBlobs, alpha, beta, d - 1, _other_player);

            #pragma omp critical
            {
                maxEval = std::max(maxEval, ev);

                alpha = std::max(alpha, ev);
                if (alpha >= beta) {
                    flag = true;
                }
            }
        }

        return maxEval;
    } else {
        Sint32 minEval = +INFINITY_BL;

        volatile bool flag = false;
        #pragma omp parallel for shared(flag)
        for (MoveEv &mv : valid_moves) {
            if (flag) continue;

            bidiarray<Sint16> currBlobs = bidiarray<Sint16>();
            copyBiDiArray(blobs, currBlobs);
            applyMove0(currBlobs, _other_player, _current_player, mv.move);
            Sint32 ev = algoAlphaBetaSortParRec(currBlobs, alpha, beta, d - 1, _current_player);

            #pragma omp critical
            {
                minEval = std::min(minEval, ev);

                beta = std::min(beta, ev);
                if (alpha >= beta) {
                    flag = true;
                }
            }
        }

        return minEval;
    }
}

movement Strategy::algoAlphaBetaSortPar(bidiarray<Sint16> &blobs, int d) {
    vector<MoveEv> valid_moves = vector<MoveEv>();
    computeValidMovesEv(blobs, _current_player, _other_player, valid_moves);
    std::sort(valid_moves.begin(), valid_moves.end(), 
               [](const MoveEv &m1, const MoveEv &m2) { return m1.ev > m2.ev; } );

    Sint32 bestEval = -INFINITY_BL;
    size_t minI = 0;
    Sint32 alpha = -INFINITY_BL;
    Sint32 beta = +INFINITY_BL;
    
    #pragma omp parallel for
    for (size_t i = 0; i < valid_moves.size(); i++) {
        bidiarray<Sint16> currBlobs = bidiarray<Sint16>();
        copyBiDiArray(blobs, currBlobs);
        applyMove0(currBlobs, _current_player, _other_player, valid_moves[i].move);
        Sint32 ev = algoAlphaBetaSortParRec(currBlobs, alpha, beta, d - 1, _other_player);

        #pragma omp critical
        {
            if (ev > bestEval || (ev == bestEval && i < minI)) {
            // if (ev > bestEval) {
                bestEval = ev;
                minI = i;
            }
        }
    }

    return valid_moves[minI].move;
}

void Strategy::algoAlphaBetaSortParAnytime() {
    int d = BASE_DEPTH;
    while (d < 100) {
        if (d <= DEPTH_LIM_0) {
            currLimDepth = d - DEPTH_VAL_0;
        } else {
            currLimDepth = d - DEPTH_VAL_1;
        }
        movement bestMove = algoAlphaBetaSortPar(_blobs, d);
        _saveBestMove(bestMove);
#ifdef DEBUG_MODE
        std::cout << "AlphaBetaSortPar ================= Depth = " << d << " Done =================" << std::endl;
#endif
        d++;
    }
}

Sint32 Strategy::algoAlphaBetaSortZobristSeqRec(bidiarray<Sint16> &blobs, Sint32 alpha, Sint32 beta, int d, Uint16 player) {
    if (d == 0) return estimateCurrentScore2(blobs);
    
    vector<movement> valid_moves = vector<movement>();
    computeValidMoves2(blobs, player, valid_moves);

    if (valid_moves.empty()) return estimateCurrentScore2(blobs);

    vector<point> changedBlobs = vector<point>();
    if (player == _current_player) {
        Sint32 maxEval = -INFINITY_BL;

        for (movement &mv : valid_moves) {
            applyMove2(blobs, _current_player, _other_player, mv, changedBlobs);
            Sint32 ev = checkZobristCallSeq(blobs, alpha, beta, d, _current_player, _other_player);
            revertMove2(blobs, _current_player, _other_player, mv, changedBlobs);

            maxEval = std::max(maxEval, ev);

            alpha = std::max(alpha, ev);
            if (alpha >= beta) break;
        }

        return maxEval;
    } else {
        Sint32 minEval = +INFINITY_BL;

        for (movement &mv : valid_moves) {
            applyMove2(blobs, _other_player, _current_player, mv, changedBlobs);
            Sint32 ev = checkZobristCallSeq(blobs, alpha, beta, d, _other_player, _current_player);
            revertMove2(blobs, _other_player, _current_player, mv, changedBlobs);

            minEval = std::min(minEval, ev);

            beta = std::min(beta, ev);
            if (alpha >= beta) break;
        }

        return minEval;
    }
}

Sint32 Strategy::algoAlphaBetaSortZobristParRec(bidiarray<Sint16> &blobs, Sint32 alpha, Sint32 beta, int d, Uint16 player) {
    if (d == 0) return estimateCurrentScore2(blobs);

    if (d <= currLimDepth) {
        return algoAlphaBetaSortZobristSeqRec(blobs, alpha, beta, d, player);
    }

    vector<MoveEv> valid_moves = vector<MoveEv>();
    computeValidMovesEv(blobs, player, (player == 1) ? 0 : 1, valid_moves);
    if (player == _current_player) {
        std::sort(valid_moves.begin(), valid_moves.end(), 
                [](const MoveEv &m1, const MoveEv &m2) { return m1.ev > m2.ev; });
    } else {
        std::sort(valid_moves.begin(), valid_moves.end(), 
                [](const MoveEv &m1, const MoveEv &m2) { return m1.ev < m2.ev; });
    }

    if (valid_moves.empty()) return estimateCurrentScore2(blobs);

    if (player == _current_player) {
        Sint32 maxEval = -INFINITY_BL;

        volatile bool flag = false;
        #pragma omp parallel for shared(flag)
        for (MoveEv &mv : valid_moves) {
            if (flag) continue;

            bidiarray<Sint16> currBlobs = bidiarray<Sint16>();
            copyBiDiArray(blobs, currBlobs);
            applyMove0(currBlobs, _current_player, _other_player, mv.move);
            Sint32 ev = checkZobristCallPar(currBlobs, alpha, beta, d, _current_player, _other_player);

            #pragma omp critical
            {
                maxEval = std::max(maxEval, ev);

                alpha = std::max(alpha, ev);
                if (alpha >= beta) {
                    flag = true;
                }
            }
        }

        return maxEval;
    } else {
        Sint32 minEval = +INFINITY_BL;

        volatile bool flag = false;
        #pragma omp parallel for shared(flag)
        for (MoveEv &mv : valid_moves) {
            if (flag) continue;

            bidiarray<Sint16> currBlobs = bidiarray<Sint16>();
            copyBiDiArray(blobs, currBlobs);
            applyMove0(currBlobs, _other_player, _current_player, mv.move);
            Sint32 ev = checkZobristCallPar(currBlobs, alpha, beta, d, _other_player, _current_player);

            #pragma omp critical
            {
                minEval = std::min(minEval, ev);

                beta = std::min(beta, ev);
                if (alpha >= beta) {
                    flag = true;
                }
            }
        }

        return minEval;
    }
}

movement Strategy::algoAlphaBetaSortZobristPar(bidiarray<Sint16> &blobs, int d) {
    vector<MoveEv> valid_moves = vector<MoveEv>();
    computeValidMovesEv(blobs, _current_player, _other_player, valid_moves);
    std::sort(valid_moves.begin(), valid_moves.end(), 
               [](const MoveEv &m1, const MoveEv &m2) { return m1.ev > m2.ev; } );

    Sint32 bestEval = -INFINITY_BL;
    size_t minI = 0;
    Sint32 alpha = -INFINITY_BL;
    Sint32 beta = +INFINITY_BL;
    
    #pragma omp parallel for
    for (size_t i = 0; i < valid_moves.size(); i++) {
        bidiarray<Sint16> currBlobs = bidiarray<Sint16>();
        copyBiDiArray(blobs, currBlobs);
        applyMove0(currBlobs, _current_player, _other_player, valid_moves[i].move);
        Sint32 ev = checkZobristCallPar(currBlobs, alpha, beta, d, _current_player, _other_player);

        #pragma omp critical
        {
            if (ev > bestEval || (ev == bestEval && i < minI)) {
            // if (ev > bestEval) {
                bestEval = ev;
                minI = i;
            }
        }
    }

    return valid_moves[minI].move;
}

void Strategy::algoAlphaBetaSortZobristParAnytime() {
    initZobrist();

    int d = 1;
    currLimDepth = 4;
    while (d < ZOBRIST_MAX_DEPTH) {
        // if (d <= DEPTH_LIM_0) {
            // currLimDepth = d - DEPTH_VAL_0;
        // } else {
            // currLimDepth = d - DEPTH_VAL_1;
        // }
        movement bestMove = algoAlphaBetaSortZobristPar(_blobs, d);
        _saveBestMove(bestMove);
#ifdef DEBUG_MODE
        std::cout << "AlphaBetaSortZobristPar ================= Depth = " << d << " Done =================" << std::endl;
#endif
        d++;
    }
}

void Strategy::initZobrist() {
    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::uniform_int_distribution<uint64_t> dist;

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            for (int k = 0; k < PIECE_COUNT; k++) {
                for (int d = 0; d < ZOBRIST_MAX_DEPTH; d++) {
                    zobristTable[i][j][k][d] = dist(gen);
                }
            }
        }
    }

    player0ToMove = dist(gen);
}

uint64_t Strategy::hashZobrist(bidiarray<Sint16> &blobs, int d, Uint16 whichPlayerTurn) {
    uint64_t h = 0;
    
    if (whichPlayerTurn == 0) {
        h = h ^ player0ToMove;
    }

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (blobs.get(i, j) != -1) {
                h = h ^ zobristTable[i][j][blobs.get(i, j)][d];
            }
        }
    }

    return h;
}

Sint32 Strategy::checkZobristCallSeq(bidiarray<Sint16> &blobs, Sint32 alpha, Sint32 beta, int d, Uint16 currPlayer, Uint16 otherPlayer) {
    uint64_t h = hashZobrist(blobs, d, currPlayer);
    if (blobsEvMap.find(h) != blobsEvMap.end()) {
// #ifdef DEBUG_MODE
        // #pragma omp critical
        // std::cout << "==================== Found ==================== " << h << std::endl;
// #endif
        return blobsEvMap[h];
    } else {
// #ifdef DEBUG_MODE
        // #pragma omp critical
        // std::cout << "==================== Not Found ==================== " << h << std::endl;
// #endif
        Sint32 ev = algoAlphaBetaSortZobristSeqRec(blobs, alpha, beta, d - 1, otherPlayer);
        #pragma omp critical
        blobsEvMap[h] = ev;
        return ev;
    }
}

Sint32 Strategy::checkZobristCallPar(bidiarray<Sint16> &blobs, Sint32 alpha, Sint32 beta, int d, Uint16 currPlayer, Uint16 otherPlayer) {
    uint64_t h = hashZobrist(blobs, d, currPlayer);
    if (blobsEvMap.find(h) != blobsEvMap.end()) {
// #ifdef DEBUG_MODE
        // #pragma omp critical
        // std::cout << "==================== Found ==================== " << h << std::endl;
// #endif
        return blobsEvMap[h];
    } else {
// #ifdef DEBUG_MODE
        // #pragma omp critical
        // std::cout << "==================== Not Found ==================== " << h << std::endl;
// #endif
        Sint32 ev = algoAlphaBetaSortZobristParRec(blobs, alpha, beta, d - 1, otherPlayer);
        #pragma omp critical
        blobsEvMap[h] = ev;
        return ev;
    }
}

Sint32 Strategy::algoPVSRec(Sint32 alpha, Sint32 beta, int d, Uint16 player) {
    if (d == 0) {
        return estimateCurrentScorePVS(player);
    }

    vector<movement> valid_moves = vector<movement>();
    computeValidMoves(player, valid_moves);

    if (valid_moves.empty()) {
        return estimateCurrentScorePVS(player);
    }

    Uint16 rival = (player == 1) ? 0 : 1;

    Sint32 ev = 0;

    vector<point> changedBlobs = vector<point>();

    for (size_t i = 0; i < valid_moves.size(); i++) {
        applyMove(player, rival, valid_moves[i], changedBlobs);
        if (i == 0) {
            ev = -algoPVSRec(-beta, -alpha, d - 1, rival);
        } else {
            ev = -algoPVSRec(-alpha - 1, -alpha, d - 1, rival);
            if (alpha < ev && ev < beta) {
                ev = -algoPVSRec(-beta, -alpha, d - 1, rival);
            }
        }
        revertMove(player, rival, valid_moves[i], changedBlobs);

        alpha = std::max(alpha, ev);
        if (alpha >= beta) {
            break;
        }
    }

    return alpha;
}

movement Strategy::algoPVS(int d) {
    vector<movement> valid_moves = vector<movement>();
    computeValidMoves(_current_player, valid_moves);

    Sint32 bestEval = -INFINITY_BL;
    size_t minI = 0;
    Sint32 alpha = -INFINITY_BL;
    Sint32 beta = +INFINITY_BL;
    Sint32 ev = 0;

    vector<point> changedBlobs = vector<point>();

    for (size_t i = 0; i < valid_moves.size(); i++) {
        applyMove(_current_player, _other_player, valid_moves[i], changedBlobs);
        if (i == 0) {
            ev = -algoPVSRec(-beta, -alpha, d - 1, _other_player);
        } else {
            ev = -algoPVSRec(-alpha - 1, -alpha, d - 1, _other_player);
            if (alpha < ev && ev < beta) {
                ev = -algoPVSRec(-beta, -alpha, d - 1, _other_player);
            }
        }
        revertMove(_current_player, _other_player, valid_moves[i], changedBlobs);

        if (ev > bestEval) {
            bestEval = ev;
            minI = i;
        }
    }

    return valid_moves[minI];
}

void Strategy::algoPVSAnytime() {
    int d = 1;
    while (d < 100) {
        movement bestMove = algoPVS(d);
        _saveBestMove(bestMove);
#ifdef DEBUG_MODE
        std::cout << "algoPVS ================= Depth = " << d << " Done =================" << std::endl;
#endif
        d++;
    }
}

Sint32 Strategy::algoPVSSortRec(Sint32 alpha, Sint32 beta, int d, Uint16 player) {
    if (d == 0) {
        return estimateCurrentScorePVS(player);
    }

    vector<movement> valid_moves = vector<movement>();
    computeValidMoves(player, valid_moves);

    if (valid_moves.empty()) {
        return estimateCurrentScorePVS(player);
    }

    sortByKillerMoves(valid_moves, currMaxDepth - d);

    Uint16 rival = (player == 1) ? 0 : 1;

    Sint32 ev = 0;

    vector<point> changedBlobs = vector<point>();

    for (size_t i = 0; i < valid_moves.size(); i++) {
        applyMove(player, rival, valid_moves[i], changedBlobs);
        if (i == 0) {
            ev = -algoPVSSortRec(-beta, -alpha, d - 1, rival);
        } else {
            ev = -algoPVSSortRec(-alpha - 1, -alpha, d - 1, rival);
            if (alpha < ev && ev < beta) {
                ev = -algoPVSSortRec(-beta, -alpha, d - 1, rival);
            }
        }
        revertMove(player, rival, valid_moves[i], changedBlobs);

        alpha = std::max(alpha, ev);
        if (alpha >= beta) {
            killerMovesSets[currMaxDepth - d].insert(valid_moves[i]);
            break;
        }
    }

    return alpha;
}

movement Strategy::algoPVSSort(int d) {
    vector<movement> valid_moves = vector<movement>();
    computeValidMoves(_current_player, valid_moves);
    sortByKillerMoves(valid_moves, 0); // Killer Moves of Move 0

    Sint32 bestEval = -INFINITY_BL;
    size_t minI = -1;
    Sint32 alpha = -INFINITY_BL;
    Sint32 beta = +INFINITY_BL;
    Sint32 ev = 0;

    vector<point> changedBlobs = vector<point>();

    for (size_t i = 0; i < valid_moves.size(); i++) {
        applyMove(_current_player, _other_player, valid_moves[i], changedBlobs);
        if (i == 0) {
            ev = -algoPVSSortRec(-beta, -alpha, d - 1, _other_player);
        } else {
            ev = -algoPVSSortRec(-alpha - 1, -alpha, d - 1, _other_player);
            if (alpha < ev && ev < beta) {
                ev = -algoPVSSortRec(-beta, -alpha, d - 1, _other_player);
            }
        }
        revertMove(_current_player, _other_player, valid_moves[i], changedBlobs);

        if (ev > bestEval) {
            bestEval = ev;
            minI = i;
        }
    }

    killerMovesSets[0].insert(valid_moves[minI]); // Move 0

    return valid_moves[minI];
}

void Strategy::algoPVSSortAnytime() {
    int d = 1;
    while (d < PVS_MAX_DEPTH) {
        currMaxDepth = d;
        movement bestMove = algoPVSSort(d);
        _saveBestMove(bestMove);
#ifdef DEBUG_MODE
        std::cout << "algoPVSSort ================= Depth = " << d << " Done =================" << std::endl;
#endif
        d++;
    }
}

void Strategy::sortByKillerMoves(vector<movement> &valid_moves, int k) {
    std::unordered_set<movement> &killerMoves = killerMovesSets[k];
    if (killerMoves.empty()) return;

    size_t n = valid_moves.size();
    int jMax = n / 6;
    int j = 0;
    for (size_t i = 0; i < n; i++) {
        if (killerMoves.find(valid_moves[i]) != killerMoves.end()) {
            movement killerMove = valid_moves[i];
            valid_moves[i] = valid_moves[j];
            valid_moves[j] = killerMove;
            j++;
            if (j == jMax) return;
        }
    }
}

Sint32 Strategy::algoPVSSortSeqRec(bidiarray<Sint16> &blobs, Sint32 alpha, Sint32 beta, int d, Uint16 player) {
    if (d == 0) return estimateCurrentScorePVS2(blobs, player);

    vector<movement> valid_moves = vector<movement>();
    computeValidMoves2(blobs, player, valid_moves);

    if (valid_moves.empty()) return estimateCurrentScorePVS2(blobs, player);

    // size_t j = sortByKillerMoves2(valid_moves, currMaxDepth - d);

    Uint16 rival = (player == 1) ? 0 : 1;

    Sint32 ev = 0;

    vector<point> changedBlobs = vector<point>();

    for (size_t i = 0; i < valid_moves.size(); i++) {
        applyMove2(blobs, player, rival, valid_moves[i], changedBlobs);
        if (i == 0) {
            ev = -algoPVSSortSeqRec(blobs, -beta, -alpha, d - 1, rival);
        } else {
            ev = -algoPVSSortSeqRec(blobs, -alpha - 1, -alpha, d - 1, rival);
            if (alpha < ev && ev < beta) {
                ev = -algoPVSSortSeqRec(blobs, -beta, -alpha, d - 1, rival);
            }
        }
        revertMove2(blobs, player, rival, valid_moves[i], changedBlobs);

        alpha = std::max(alpha, ev);
        if (alpha >= beta) {
            // #pragma omp critical
            // if (i >= j) killerMovesSets[currMaxDepth - d].insert(valid_moves[i]);
            break;
        }
    }

    return alpha;
}

Sint32 Strategy::algoPVSSortParRec(bidiarray<Sint16> &blobs, Sint32 alpha, Sint32 beta, int d, Uint16 player) {
    if (d == 0) return estimateCurrentScorePVS2(blobs, player);

    if (d <= currLimDepth) {
        return algoPVSSortSeqRec(blobs, alpha, beta, d, player);
    }

    vector<movement> valid_moves = vector<movement>();
    computeValidMoves2(blobs, player, valid_moves);

    if (valid_moves.empty()) return estimateCurrentScorePVS2(blobs, player);

    size_t j = sortByKillerMoves2(valid_moves, currMaxDepth - d);

    Uint16 rival = (player == 1) ? 0 : 1;

    volatile bool flag = false;

//////////////////////////////////////////////////////////////////////////// No Parallel Start
    vector<point> changedBlobs = vector<point>();
    for (size_t i = 0; i < j; i++) {
        Sint32 ev = 0;
        applyMove2(blobs, player, rival, valid_moves[i], changedBlobs);
        if (i == 0) {
            ev = -algoPVSSortParRec(blobs, -beta, -alpha, d - 1, rival);
        } else {
            ev = -algoPVSSortParRec(blobs, -alpha - 1, -alpha, d - 1, rival);
            if (alpha < ev && ev < beta) {
                ev = -algoPVSSortParRec(blobs, -beta, -alpha, d - 1, rival);
            }
        }
        revertMove2(blobs, player, rival, valid_moves[i], changedBlobs);

        alpha = std::max(alpha, ev);
        if (alpha >= beta) {
            goto ret;
        }
    }
////////////////////////////////////////////////////////////////////////////  No Parallel End

//////////////////////////////////////////////////////////////////////////// Parallel Start
    #pragma omp parallel for shared(flag)
    for (size_t i = j; i < valid_moves.size(); i++) {
        if (flag) continue;

        Sint32 ev = 0;
        bidiarray<Sint16> currBlobs = bidiarray<Sint16>();
        copyBiDiArray(blobs, currBlobs);
        applyMove0(currBlobs, player, rival, valid_moves[i]);
        if (i == 0) {
            ev = -algoPVSSortParRec(currBlobs, -beta, -alpha, d - 1, rival);
        } else {
            ev = -algoPVSSortParRec(currBlobs, -alpha - 1, -alpha, d - 1, rival);
            if (alpha < ev && ev < beta) {
                ev = -algoPVSSortParRec(currBlobs, -beta, -alpha, d - 1, rival);
            }
        }

        #pragma omp critical
        {
            alpha = std::max(alpha, ev);
            if (alpha >= beta) {
                killerMovesSets[currMaxDepth - d].insert(valid_moves[i]);
                flag = true;
            }
        }
    }
//////////////////////////////////////////////////////////////////////////// Parallel End

ret:
    return alpha;
}

movement Strategy::algoPVSSortPar(bidiarray<Sint16> &blobs, int d) {
    if (d <= 2) {
        return algoPVSSort(d);
    }

    vector<movement> valid_moves = vector<movement>();
    computeValidMoves2(blobs, _current_player, valid_moves);

    size_t j = sortByKillerMoves2(valid_moves, 0); // Killer Moves of Move 0

    Sint32 bestEval = -INFINITY_BL;
    size_t minI = 42424242;
    Sint32 alpha = -INFINITY_BL;
    Sint32 beta = +INFINITY_BL;

//////////////////////////////////////////////////////////////////////////// No Parallel Start
    vector<point> changedBlobs = vector<point>();
    for (size_t i = 0; i < j; i++) {
        Sint32 ev = 0;
        applyMove2(blobs, _current_player, _other_player, valid_moves[i], changedBlobs);
        if (i == 0) {
            ev = -algoPVSSortParRec(blobs, -beta, -alpha, d - 1, _other_player);
        } else {
            ev = -algoPVSSortParRec(blobs, -alpha - 1, -alpha, d - 1, _other_player);
            if (alpha < ev && ev < beta) {
                ev = -algoPVSSortParRec(blobs, -beta, -alpha, d - 1, _other_player);
            }
        }
        revertMove2(blobs, _current_player, _other_player, valid_moves[i], changedBlobs);

        if (ev > bestEval) {
            bestEval = ev;
            minI = i;
        }
    }
////////////////////////////////////////////////////////////////////////////  No Parallel End

//////////////////////////////////////////////////////////////////////////// Parallel Start
    #pragma omp parallel for
    for (size_t i = j; i < valid_moves.size(); i++) {
        Sint32 ev = 0;
        bidiarray<Sint16> currBlobs = bidiarray<Sint16>();
        copyBiDiArray(blobs, currBlobs);
        applyMove0(currBlobs, _current_player, _other_player, valid_moves[i]);
        if (i == 0) {
            ev = -algoPVSSortParRec(currBlobs, -beta, -alpha, d - 1, _other_player);
        } else {
            ev = -algoPVSSortParRec(currBlobs, -alpha - 1, -alpha, d - 1, _other_player);
            if (alpha < ev && ev < beta) {
                ev = -algoPVSSortParRec(currBlobs, -beta, -alpha, d - 1, _other_player);
            }
        }

        #pragma omp critical
        {
            if (ev > bestEval || (ev == bestEval && i < minI)) {
                bestEval = ev;
                minI = i;
            }
        }
    }
//////////////////////////////////////////////////////////////////////////// Parallel End

    killerMovesSets[0].insert(valid_moves[minI]); // Move 0

    return valid_moves[minI];
}

void Strategy::algoPVSSortParAnytime() {
    // omp_init_lock(&writelock);

    int d = 1;
    while (d < PVS_MAX_DEPTH) {
        currMaxDepth = d;
        if (d <= DEPTH_LIM_0) {
            currLimDepth = d - DEPTH_VAL_0;
        } else {
            currLimDepth = d - DEPTH_VAL_1;
        }
        movement bestMove = algoPVSSortPar(_blobs, d);
        _saveBestMove(bestMove);
#ifdef DEBUG_MODE
        std::cout << "algoPVSSortPar ================= Depth = " << d << " Done =================" << std::endl;
#endif
        d++;
    }
}

size_t Strategy::sortByKillerMoves2(vector<movement> &valid_moves, int k) {
    std::unordered_set<movement> &killerMoves = killerMovesSets[k];
    if (killerMoves.empty()) return 0;

    size_t n = valid_moves.size();
    size_t jMax = n / 6;
    size_t j = 0;
    bool v = false;
    for (size_t i = 0; i < n; i++) {
        #pragma omp critical
        v = killerMoves.find(valid_moves[i]) != killerMoves.end();
        if (v) {
            movement killerMove = valid_moves[i];
            valid_moves[i] = valid_moves[j];
            valid_moves[j] = killerMove;
            j++;
            if (j == jMax) return j;
        }
    }

    return j;
}

void Strategy::computeBestMove() {
#ifdef DEBUG_MODE
    std::cout << "Player : " << _current_player << std::endl;
#endif

#ifdef DEBUG_MODE
    auto startTime = std::chrono::high_resolution_clock::now();
#endif

#ifdef DEBUG_MODE
    if (_current_player == 1) { // Blue
        // algoNaif();
        // algoGlouton();
        // algoMinMax(4);
        // algoAlphaBeta(6);
        // algoAlphaBetaPar(_blobs, 8);
        // algoMinMax2Anytime();
        // algoAlphaBeta2Anytime();
        // algoAlphaBetaPar2Anytime();
        // algoAlphaBetaSortParAnytime();
        // algoAlphaBetaSortZobristParAnytime();
        // algoPVSAnytime();
        // algoPVSSortAnytime();
        algoPVSSortParAnytime();
    } else { // Red
        // algoNaif();
        // algoGlouton();
        // algoMinMax(4);
        // algoAlphaBeta(4);
        // algoAlphaBetaPar(_blobs, 8);
        // algoMinMax2Anytime();
        // algoAlphaBeta2Anytime();
        // algoAlphaBetaPar2Anytime();
        algoAlphaBetaSortParAnytime();
        // algoAlphaBetaSortZobristParAnytime();
        // algoPVSAnytime();
        // algoPVSSortAnytime();
        // algoPVSSortParAnytime();
    }
#else
    algoPVSSortParAnytime();
#endif

#ifdef DEBUG_MODE
    auto endTime = std::chrono::high_resolution_clock::now();
    std::chrono::milliseconds t = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);

    std::cout << "TIME : " << t.count() / 1000.0 << std::endl;
#endif
}
