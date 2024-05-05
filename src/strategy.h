#ifndef __STRATEGY_H
#define __STRATEGY_H

#include "common.h"
#include "bidiarray.h"
#include "move.h"
#include "moveEv.h"
#include <unordered_map>
#include <unordered_set>
#include <omp.h>

#define INFINITY_BL 96 // > N * N = 64
#define N 8
#define PIECE_COUNT 2
#define ZOBRIST_MAX_DEPTH 32
#define PVS_MAX_DEPTH 20

class Strategy {

private:
    //! array containing all blobs on the board
    bidiarray<Sint16> _blobs;
    //! an array of booleans indicating for each cell whether it is a hole or not.
    const bidiarray<bool> &_holes;
    //! Current player
    Uint16 _current_player;
    
    Uint16 _other_player;
    
    //! Call this function to save your best move.
    //! Multiple call can be done each turn,
    //! Only the last move saved will be used.
    void (*_saveBestMove)(movement &);

    std::unordered_map<uint64_t, Sint32> blobsEvMap;
    uint64_t zobristTable[N][N][PIECE_COUNT][ZOBRIST_MAX_DEPTH];
    uint64_t player0ToMove;

    int currLimDepth;

    int currMaxDepth;
    std::unordered_set<movement> killerMovesSets[PVS_MAX_DEPTH];

    // omp_lock_t writelock;

public:
    Strategy(bidiarray<Sint16> &blobs, 
              const bidiarray<bool> &holes,
              const Uint16 current_player,
              void (*saveBestMove)(movement &))
            : _blobs(blobs),_holes(holes), _current_player(current_player),
              _other_player((current_player == 1) ? 0 : 1), _saveBestMove(saveBestMove),
              blobsEvMap(), zobristTable(), player0ToMove(), currLimDepth(),
              currMaxDepth(), killerMovesSets() {

    }

    Strategy(const Strategy& St)
            : _blobs(St._blobs), _holes(St._holes), _current_player(St._current_player),
            _other_player((St._current_player == 1) ? 0 : 1), _saveBestMove(St._saveBestMove),
              blobsEvMap(), zobristTable(), player0ToMove(), currLimDepth(),
              currMaxDepth(), killerMovesSets() {

    }

    ~Strategy() { }

    void algoNaif();
    
    void algoGlouton();

    Sint32 algoMinMaxRec(int d, Uint16 player);

    void algoMinMax(int d);

    movement algoMinMax2(int d);

    void algoMinMax2Anytime();

    Sint32 algoAlphaBetaRec(Sint32 alpha, Sint32 beta, int d, Uint16 player);

    void algoAlphaBeta(int d);

    movement algoAlphaBeta2(int d);

    void algoAlphaBeta2Anytime();

    inline int getMoveType(movement &mv) { // 1 si une case, 0 sinon.
          return std::abs(mv.nx - mv.ox) < 2 &&
           std::abs(mv.ny - mv.oy) < 2;
    }

    void applyMove(Uint16 ally, Uint16 rival, movement &mv, vector<point> &changedBlobs);

    void revertMove(Uint16 ally, Uint16 rival, movement &mv, vector<point> &changedBlobs);

    vector<movement> &computeValidMoves(Uint16 player, vector<movement> &valid_moves) const;

    Sint32 estimateCurrentScore() const;

    Sint32 algoAlphaBetaSeqRec(bidiarray<Sint16> &blobs, Sint32 alpha, Sint32 beta, int d, Uint16 player);

    Sint32 algoAlphaBetaParRec(bidiarray<Sint16> &blobs, Sint32 alpha, Sint32 beta, int d, Uint16 player);

    void algoAlphaBetaPar(bidiarray<Sint16> &blobs, int d);

    inline void copyBiDiArray(bidiarray<Sint16> &src, bidiarray<Sint16> &dst) {
      memcpy(dst.array, src.array, N * N * sizeof(Sint16));
    }

    void applyMove0(bidiarray<Sint16> &blobs, Uint16 ally, Uint16 rival, movement &mv);

    void applyMove2(bidiarray<Sint16> &blobs, Uint16 ally, Uint16 rival, movement &mv, vector<point> &changedBlobs);

    void revertMove2(bidiarray<Sint16> &blobs, Uint16 ally, Uint16 rival, movement &mv, vector<point> &changedBlobs, bool doClear = true);

    vector<movement> &computeValidMoves2(bidiarray<Sint16> &blobs, Uint16 player, vector<movement> &valid_moves) const;

    Sint32 estimateCurrentScore2(bidiarray<Sint16> &blobs) const;

    movement algoAlphaBetaPar2(bidiarray<Sint16> &blobs, int d);

    void algoAlphaBetaPar2Anytime();

    void computeValidMovesEv(bidiarray<Sint16> &blobs, Uint16 ally, Uint16 rival, vector<MoveEv> &valid_moves);

    Sint32 algoAlphaBetaSortParRec(bidiarray<Sint16> &blobs, Sint32 alpha, Sint32 beta, int d, Uint16 player);

    movement algoAlphaBetaSortPar(bidiarray<Sint16> &blobs, int d);

    void algoAlphaBetaSortParAnytime();

    Sint32 algoAlphaBetaSortZobristSeqRec(bidiarray<Sint16> &blobs, Sint32 alpha, Sint32 beta, int d, Uint16 player);

    Sint32 algoAlphaBetaSortZobristParRec(bidiarray<Sint16> &blobs, Sint32 alpha, Sint32 beta, int d, Uint16 player);

    movement algoAlphaBetaSortZobristPar(bidiarray<Sint16> &blobs, int d);

    void algoAlphaBetaSortZobristParAnytime();

    void initZobrist();

    uint64_t hashZobrist(bidiarray<Sint16> &blobs, int d, Uint16 whichPlayerTurn);

    Sint32 checkZobristCallSeq(bidiarray<Sint16> &blobs, Sint32 alpha, Sint32 beta, int d, Uint16 currPlayer, Uint16 otherPlayer);

    Sint32 checkZobristCallPar(bidiarray<Sint16> &blobs, Sint32 alpha, Sint32 beta, int d, Uint16 currPlayer, Uint16 otherPlayer);

    Sint32 algoPVSRec(Sint32 alpha, Sint32 beta, int d, Uint16 player);

    movement algoPVS(int d);

    void algoPVSAnytime();

    inline Sint32 estimateCurrentScorePVS(Uint16 player) {
      if (player == _current_player) return estimateCurrentScore();
      else return -estimateCurrentScore();
    }

    Sint32 algoPVSSortRec(Sint32 alpha, Sint32 beta, int d, Uint16 player);

    movement algoPVSSort(int d);

    void algoPVSSortAnytime();

    void sortByKillerMoves(vector<movement> &valid_moves, int k);

    Sint32 algoPVSSortSeqRec(bidiarray<Sint16> &blobs, Sint32 alpha, Sint32 beta, int d, Uint16 player);

    Sint32 algoPVSSortParRec(bidiarray<Sint16> &blobs, Sint32 alpha, Sint32 beta, int d, Uint16 player);

    movement algoPVSSortPar(bidiarray<Sint16> &blobs, int d);

    void algoPVSSortParAnytime();

    inline Sint32 estimateCurrentScorePVS2(bidiarray<Sint16> &blobs, Uint16 player) const {
        if (player == _current_player) return estimateCurrentScore2(blobs);
        else return -estimateCurrentScore2(blobs);
    }

    size_t sortByKillerMoves2(vector<movement> &valid_moves, int k);

    void computeBestMove();
};

#endif
