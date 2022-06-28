/************************************************************
 * HMMER - Biological sequence analysis with HMMs
 * Copyright 1992-1995 Sean R. Eddy
 *
 *   This source code is distributed under the terms of the
 *   GNU General Public License. See the files COPYING and
 *   GNULICENSE for details.
 *
 ************************************************************/

/* weeviterbi.c
 * SRE, Thu Nov 10 11:13:08 1994
 *
 * Alignment of a complete sequence to all or part of an HMM
 * [modified Smith/Waterman] in linear memory. Alignment is
 * done by a divide-and-conquer method.
 * 
 * In Smith/Waterman mode, WeeViterbi() uses the same probabilistic 
 * model as SWViterbi(). It is controlled by two parameters P2 and
 * P3. P2 is the probability of entering at match state 1. P3 is the 
 * probability of leaving from any internal match state. P1, which 
 * controls the scoring of external insertions in SWViterbi(), is 
 * not used here because the full sequence must align.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <string.h>
#include <math.h>

#include "states.h"
#include "squid.h"
#include "externs.h"

#ifdef MEMDEBUG
#include "dbmalloc.h"
#endif

static int *startcost;         /* position specific scores for entering S/W match */
static int *endcost;           /* position specific scores for leaving S/W match  */

/* Function: get_midpt()
 * 
 * Purpose:  Divide and conquer alignment -- run forwards and backwards
 *           over the dynamic programming matrix, find and return an
 *           optimal midpoint. This would be called recursively on
 *           progressively smaller bits of the sequence and model to
 *           construct a complete alignment.
 *           
 * Args:     shmm   - model, integer log odds form
 *           seq    - 1..L prepared sequence      
 *           L      - length of sequence
 *           smithwaterman - TRUE if we're allowing internal starts/ends on model
 *           k1     - model node to start at, 0..M+1
 *           t1     - state type to start at, MATCH or INSERT
 *           s1     - sequence position to start at, 0..L+1
 *           k3     - model node to end at, 0..M+1
 *           t3     - state type to end at, MATCH or INSERT
 *           s3     - sequence position to end at, 0..L+1
 *           ret_k2 - RETURN: midpoint, node position in model
 *           ret_t2 - RETURN: midpoint, state type, MATCH or INSERT
 *           ret_s2 - RETURN: midpoint, sequence position, 0..L+1
 *           ret_max- RETURN: best score at the midpoint
 */
static void 
get_midpt(struct shmm_s *shmm, char *seq, int L, int smithwaterman,
	  int k1, int t1, int s1, int k3, int t3, int s3,
	  int *ret_k2, int *ret_t2, int *ret_s2, int *ret_max)
{
  struct vit_s *mx[2];          /* two scoring rows, for forwards/backwards */
  struct vit_s *fwd;            /* saved final forward row at midpt         */
  struct vit_s *cur;            /* pointer to current row of mx             */
  struct vit_s *nxt;            /* pointer to next row of mx                */
  struct vit_s *init_row;       /* pre-stored row for erasing               */
  int init_row_size;		/* length of init_row in bytes              */
  int i;			/* counter for sequence position: 0,1..L    */
  int k;			/* counter for model node position: 0,1..M  */
  int score;			/* tmp for holding intermediate score       */
  int k2, t2, s2;		/* RETURN: info on alignment midpoint       */
  int i_symidx;			/* index of sequence character              */
  int *m_emit;                  /* MAT state emission score vector          */
  int *i_emit;                  /* INS state emission score vector          */
  int *tptr;                    /* state transition score vector            */
  int  max;			/* best score in midpoint row               */
  
  /* Choose a midpoint on the sequence, s2
   */
  s2 = s1 + (s3-s1) / 2;
  
  /********************************************
   * Initialization
   ********************************************/
				/* allocate a scrolling window for the calculation matrix,
				   which is 2 rows by 0..M+2 cols */
  for (i = 0; i < 2 ; i++)
    if ((mx[i] = (struct vit_s *) malloc ((shmm->M+3) * sizeof(struct vit_s))) == NULL) 
      Die("malloc failed"); 
 
  /* create init_row, an initialized row of scores
   * that can be slapped in place by a memcpy() 
   * need a dummy column to the right, as well as the "dummy" rows
   * 0 and M+1, to deal with boundary conditions of back pass
   */
  init_row_size = (shmm->M + 3) * sizeof(struct vit_s);
  if ((init_row = (struct vit_s *) malloc (init_row_size)) == NULL ||
      (fwd      = (struct vit_s *) malloc (init_row_size)) == NULL)
    Die("malloc failed");
  for (k = 0; k <= shmm->M + 1; k++)
    {
      init_row[k].score_m = -99999999;
      init_row[k].score_d = -99999999;
      init_row[k].score_i = -99999999;
    }

  /********************************************
   * Forward pass from start to midpoint on sequence
   ********************************************/
				/* initialize the first row */
  memcpy(mx[s1%2], init_row, init_row_size);
  switch (t1) {
  case MATCH:  mx[s1%2][k1].score_m = 0; break;
  case INSERT: mx[s1%2][k1].score_i = 0; break;
  }
  for (i = s1; i < s2; i++)
    {
      cur = mx[i%2];
      nxt = mx[(i+1)%2];
				/* erase next row, set scoring ptrs */
      memcpy(nxt, init_row, init_row_size);
      i_symidx = seq[i] - 'A';
      m_emit = shmm->m_emit[i_symidx] + k1;
      i_emit = shmm->i_emit[i_symidx] + k1;
      tptr   = shmm->t + 9 * k1;

      for (k = k1; k <= k3 && k != shmm->M+1; k++)
	{
				/* add emission scores to current cell */
	  if (i > 0)
	    {
	      cur[k].score_m += *m_emit; m_emit++;
	      cur[k].score_i += *i_emit; i_emit++;
	    }
				/* deal with transitions out of delete state */
				/* note: tptr assumes order dd,di,dm,id,ii,im,md,mi,mm */
				/* we initialize with these */
	  cur[k+1].score_d = cur[k].score_d + *tptr; tptr++;
	  nxt[k].score_i   = cur[k].score_d + *tptr; tptr++;
	  nxt[k+1].score_m = cur[k].score_d + *tptr; tptr++;

				/* deal with transitions out of insert state */
	  if ((score = cur[k].score_i + *tptr) > cur[k+1].score_d)
	    cur[k+1].score_d = score;
	  tptr++;
	  if ((score = cur[k].score_i + *tptr) > nxt[k].score_i)
	    nxt[k].score_i = score;
	  tptr++;
	  if ((score = cur[k].score_i + *tptr) > nxt[k+1].score_m)
	    nxt[k+1].score_m = score;
	  tptr++;

				/* deal with transitions out of match state */
	  if ((score = cur[k].score_m + *tptr) > cur[k+1].score_d)
	    cur[k+1].score_d = score;
	  tptr++;
	  if ((score = cur[k].score_m + *tptr) > nxt[k].score_i)
	    nxt[k].score_i = score;
	  tptr++;
	  if ((score = cur[k].score_m + *tptr) > nxt[k+1].score_m)
	    nxt[k+1].score_m = score;
	  tptr++;
				/* deal with Smith/Waterman starts */
	  if (smithwaterman && i == 0)
	    if (startcost[k+1] > nxt[k+1].score_m)
	      nxt[k+1].score_m = startcost[k+1];
	}
    }
  /* Now, we have row s2 (the midpoint) as nxt. Its DEL scores are
   * not valid, but its MAT and INS scores are, and that's all we need.
   * Save it for later.
   */
  memcpy(fwd, nxt, init_row_size);

  /********************************************
   * Backward pass from end to midpoint on sequence
   * Back pass counts emission score at midpoint;
   * forward pass didn't
   ********************************************/
				/* initialize the last row */
  memcpy(mx[s3%2], init_row, init_row_size);
  switch (t3) {
  case MATCH:  mx[s3%2][k3].score_m = 0; break;
  case INSERT: mx[s3%2][k3].score_i = 0; break;
  }

  for (i = s3-1; i >= s2; i--)
    {
      cur = mx[i%2];
      nxt = mx[(i+1)%2];	
				/* erase current row, set scoring ptrs */
      memcpy(cur, init_row, init_row_size);
      i_symidx = seq[i] - 'A';
      m_emit   = shmm->m_emit[i_symidx] + k3; 
      i_emit   = shmm->i_emit[i_symidx] + k3;
      tptr     = shmm->t + ((k3+1) * 9) - 1;

      for (k = k3; k >= k1; k--)
	if (k <= shmm->M)	/* *tptr invalid for M+1! */
	  {
				/* Find best way into current match state */
				/* note: tptr assumes order dd,di,dm,id,ii,im,md,mi,mm */
	    cur[k].score_m = nxt[k+1].score_m + *tptr; 
	    tptr--; 
	    if ((score = nxt[k].score_i + *tptr) > cur[k].score_m)
	      cur[k].score_m = score;
	    tptr--;
	    if ((score = cur[k+1].score_d + *tptr) > cur[k].score_m)
	      cur[k].score_m = score;
	    tptr--;
				/* Find best way into current insert state */
	    cur[k].score_i = nxt[k+1].score_m + *tptr; 
	    tptr--;
	    if ((score = nxt[k].score_i + *tptr) > cur[k].score_i)
	      cur[k].score_i = score;
	    tptr--;
	    if ((score = cur[k+1].score_d + *tptr) > cur[k].score_i)
	      cur[k].score_i = score;
	    tptr--;
				/* Find best way into current delete state */
	    cur[k].score_d = nxt[k+1].score_m + *tptr; 
	    tptr--;
	    if ((score = nxt[k].score_i + *tptr) > cur[k].score_d)
	      cur[k].score_d = score;
	    tptr--;
	    if ((score = cur[k+1].score_d + *tptr) > cur[k].score_d)
	      cur[k].score_d = score;
	    tptr--;
				/* deal with Smith/Waterman ends */
	    if (smithwaterman && i == L)
	      if (endcost[k] > cur[k].score_m)
		cur[k].score_m = endcost[k];
	    
  			/* add in emission scores to the current cell. */
	    cur[k].score_m += *m_emit; m_emit--;
	    cur[k].score_i += *i_emit; i_emit--;
	  }
	else			/* state M+1 */
	  {
	    tptr -= 9;
	    m_emit--;
	    i_emit--;
	  }
    }

  /* cur now holds the appropriate backward-calculated row for the
   * midpoint at s2. Sum across this row, choose best midpoint.
   */
  max = -99999999;
  for (k = k1; k <= k3; k++)
    {
      if ((score = fwd[k].score_m + cur[k].score_m) > max)
	{ k2 = k; t2 = MATCH; max = score; }
      if ((score = fwd[k].score_i + cur[k].score_i) > max)
	{ k2 = k; t2 = INSERT; max = score; }
    }

  /********************************************
   * Garbage collection and return
   ********************************************/
  free(mx[0]);
  free(mx[1]);
  free(init_row);
  free(fwd);

  *ret_k2 = k2;
  *ret_t2 = t2;
  *ret_s2 = s2;
  *ret_max= max;
  return;
}



/* Function: WeeViterbi()
 * 
 * Purpose:  Alignment of an HMM to a sequence in linear memory.
 *           Can do either global, Needleman-Wunsch style alignment,
 *           or (if smithwaterman is TRUE), can do full sequence/
 *           partial model alignments. This is Smith/Waterman in
 *           a very restricted sense -- the whole sequence must match!
 *           This is used, for instance, when alignments are
 *           reconstructed for the linear-memory scanning Smith-Waterman 
 *           functions FragViterbi() and WeeSWViterbi().
 *           
 * Args:     shmm          - model, in log-odds search form
 *           seq           - prepared sequence 1..L
 *           L             - length of seq
 *           smithwaterman - TRUE to allow fragment of model to match seq
 *           P2            - S/W parameter controlling entry into model
 *           P3            - S/W parameter controlling exit from model
 *           ret_tr        - RETURN: traceback of alignment
 *           ret_score     - RETURN: score of alignment
 */
void
WeeViterbi(struct shmm_s *shmm, char *seq, int L, 
	   int smithwaterman, float P2, float P3,
	   struct trace_s **ret_tr, float *ret_sc)
{
  struct trace_s *tr;           /* RETURN: traceback */
  int *kassign;                 /* 0..L+1, alignment of seq positions to model nodes */
  int *tassign;                 /* 0..L+1, alignment of seq positions to state types */
  int *endlist;                 /* stack of end points on sequence */
  int *startlist;               /* stack of start points on sequence */
  int  lpos;			/* position in lists */
  int k2, t2, s2;		/* midpoint */
  int k1, t1, s1;		/* start point */
  int k3, t3, s3;		/* end point  */
  int score;
  int delcount, del;

  /* Set up Smith/Waterman start/stop scores, if necessary
   */
  if (smithwaterman)
    {
      startcost = (int *) MallocOrDie (sizeof(int) * (shmm->M+2));
      endcost   = (int *) MallocOrDie (sizeof(int) * (shmm->M+2));

      for (k2 = 2; k2 <= shmm->M; k2++)
	startcost[k2] = (int) (INTSCALE * LOG2((1.0-P2)/(double)(shmm->M-1)));
      startcost[0]         = -99999999;
      startcost[1]         = (int) (INTSCALE * LOG2(P2));
      startcost[shmm->M+1] = -99999999;
      
      for (k2 = 1; k2 < shmm->M; k2++)
	endcost[k2] = (int) (INTSCALE * LOG2(P3));
      endcost[0]         = -99999999;
      endcost[shmm->M]   = 0;
      endcost[shmm->M+1] = -99999999;
    }

  /* Initialize for recursive alignment
   */
  kassign   = (int *) MallocOrDie (sizeof(int) * (L+2));
  tassign   = (int *) MallocOrDie (sizeof(int) * (L+2));
  endlist   = (int *) MallocOrDie (sizeof(int) * (L+2));
  startlist = (int *) MallocOrDie (sizeof(int) * (L+2));
  
  kassign[0]   = 0;
  kassign[L+1] = shmm->M+1;
  tassign[0] = tassign[L+1] = MATCH;
  lpos      = 0;
  startlist[lpos] = 0;
  endlist[lpos]   = L+1;

  /* Recursive divide-and-conquer alignment
   */
  while (lpos >= 0)
    {
      s1 = startlist[lpos];
      k1 = kassign[s1];
      t1 = tassign[s1];
      s3 = endlist[lpos];
      k3 = kassign[s3];
      t3 = tassign[s3];
      lpos--;

      get_midpt(shmm, seq, L, smithwaterman, k1, t1, s1, k3, t3, s3, &k2, &t2, &s2, &score);
      kassign[s2] = k2;
      tassign[s2] = t2;
				/* score is valid on first pass */
      if (s1 == 0 && s3 == L+1) *ret_sc = ((float) score) / INTSCALE;

      if (s2 - s1 > 1)
	{
	  lpos++;
	  startlist[lpos] = s1;
	  endlist[lpos]   = s2;
	}
      if (s3 - s2 > 1)
	{
	  lpos++;
	  startlist[lpos] = s2;
	  endlist[lpos]   = s3;
	}
    }

  /* Construct a traceback structure
   */
  AllocTrace(L + shmm->M + 3, &tr);
  tr->tlen = 0;
  
  for (s2 = 0; s2 <= L+1; s2++)
    {
      tr->nodeidx[tr->tlen]   = kassign[s2];
      tr->statetype[tr->tlen] = tassign[s2];
      tr->rpos[tr->tlen]      = s2 -1;          /* convert to 0..L-1 */
      if (s2 == L+1) tr->rpos[tr->tlen] = -1;
      tr->tlen++;

      /* Add delete states to traceback, if necessary 
       * Smith/Waterman must be careful not to add deletes at start/end
       */
      if (s2 != L+1 &&
	  ! (smithwaterman && s2 == 0) &&
	  ! (smithwaterman && s2 == L))
	{
	  delcount = kassign[s2+1] - kassign[s2] - 1;
	  if (tassign[s2+1] == INSERT) delcount++;

	  k2 = kassign[s2] + 1;
	  for (del = 0; del < delcount; del++)
	    {
	      tr->nodeidx[tr->tlen]   = k2;
	      tr->statetype[tr->tlen] = DELETE;
	      tr->rpos[tr->tlen]      = -1;
	      tr->tlen++;
	      k2++;
	    }
	}
    }

  if (smithwaterman)
    {
      free(startcost);
      free(endcost);
    }
  free(kassign);
  free(tassign);
  free(startlist);
  free(endlist);
  *ret_tr = tr;
}

