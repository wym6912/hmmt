/************************************************************
 * HMMER - Biological sequence analysis with HMMs
 * Copyright 1992-1995 Sean R. Eddy
 *
 *   This source code is distributed under the terms of the
 *   GNU General Public License. See the files COPYING and
 *   GNULICENSE for details.
 *
 ************************************************************/

/* viterbi.c
 * optimized according to ideas from James Crook, Mon Mar  1 11:30:51 1993
 * more heavy optimization Tue Jul  5 16:04:33 1994
 *
 * Implementation of the Viterbi dynamic programming
 * algorithm for aligning an HMM (of length M) to a sequence (of
 * length L).
 * 
 * ViterbiFill() fills in an 0..L+1 by 0..M+1 matrix using the Viterbi
 * dynamic programming algorithm, and returns 1 on success,
 * 0 on failure. The best global alignment score is mx[L+1][M+1].score_m.
 *
 * The best score is the maximum score.
 *
 * The caller may only be interested in the score,
 * but it must remember to free the memory in mx.
 * 
 * ViterbiTrace() then takes that matrix and traces back through
 * it, to find the optimal state sequence. The state sequence
 * is returned in a traceback structure trace_s. tr->nodeidx stores the index (k)
 * and tr->statetype stores the kind of state (insert, delete,
 * or match). Be careful here... both of these are arrays of 0..N-1,
 * where N is the *number of states in the optimal sequence*,
 * not the number of states in the model! Also remember that the
 * optimal state sequence will always include two dummy states,
 * for BEGIN and END, at position 0 and N-1, with nodeidx's of
 * 0 and M+1;  a la Haussler, both are represented by MATCH states.
 * This is marginally unsatisfactory, as we must make sure
 * they never get penalized for not emitting a symbol.
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#ifdef DEBUG
#include <assert.h>
#endif 

#include "states.h"
#include "externs.h"
#include "squid.h"

#ifdef MEMDEBUG
#include "dbmalloc.h"
#endif


/* Function: MakeSearchHMM()
 * 
 * Purpose:  Convert an HMM (probability form) to integer log-odds 
 *           form, for searching/alignment algorithms.
 *           
 * Args:     hmm       - probability-form hmm to convert
 *           randomseq - random sequence model
 *           shmm      - integer log-odds search HMM to create (pre-allocated)
 *
 * Return:   (void)
 */
void
MakeSearchHMM(struct hmm_struc *hmm, float *randomseq, struct shmm_s *shmm)
{
  int              k, sym, x;
  float            tmp;

  /* Symbol emission probabilities.
   * The search model keeps vectors for all 26 letters, for speedy lookup.
   */
  for (sym = 'A'; sym <= 'Z'; sym++)
    {
      if (strchr(Alphabet, sym) != NULL)
	{
	  x = SYMIDX(sym);
	  for (k = 0; k <= hmm->M; k++)
	    {
	      tmp = LOG2(hmm->mat[k].p[x] / randomseq[x]);
	      shmm->m_emit[sym-'A'][k] = (int) (INTSCALE * tmp);

				/* Inserts are HARDWIRED to give zero score!!! */
	      shmm->i_emit[sym-'A'][k] = 0;
#ifdef SRE_REMOVED
	      tmp = LOG2(hmm->ins[k].p[x] / prior->p[INSERT][x]);
	      shmm->i_emit[sym-'A'][k] = (int) (INTSCALE * tmp);
#endif
	    }
	}
      else		/* degenerate symbols */
	{
	  for (k = 0; k <= hmm->M; k++)
	    {
	      shmm->m_emit[sym-'A'][k] = Symscore(sym, hmm->mat[k].p, randomseq);
	      shmm->i_emit[sym-'A'][k] = Symscore(sym, hmm->ins[k].p, randomseq);
	    }
	}
    }

  /* State transition probabilities
   */
  for (k = 0; k <= hmm->M; k++)
    {
      tmp = LOG2(hmm->del[k].t[DELETE]); shmm->t[k*9 + Tdd] = (int) (INTSCALE * tmp); 
      tmp = LOG2(hmm->del[k].t[INSERT]); shmm->t[k*9 + Tdi] = (int) (INTSCALE * tmp); 
      tmp = LOG2(hmm->del[k].t[MATCH]);  shmm->t[k*9 + Tdm] = (int) (INTSCALE * tmp); 

      tmp = LOG2(hmm->ins[k].t[DELETE]); shmm->t[k*9 + Tid] = (int) (INTSCALE * tmp); 
      tmp = LOG2(hmm->ins[k].t[INSERT]); shmm->t[k*9 + Tii] = (int) (INTSCALE * tmp); 
      tmp = LOG2(hmm->ins[k].t[MATCH]);  shmm->t[k*9 + Tim] = (int) (INTSCALE * tmp); 

      tmp = LOG2(hmm->mat[k].t[DELETE]); shmm->t[k*9 + Tmd] = (int) (INTSCALE * tmp); 
      tmp = LOG2(hmm->mat[k].t[INSERT]); shmm->t[k*9 + Tmi] = (int) (INTSCALE * tmp); 
      tmp = LOG2(hmm->mat[k].t[MATCH]);  shmm->t[k*9 + Tmm] = (int) (INTSCALE * tmp); 
    }

  /* Annotation
   */
  shmm->flags = hmm->flags;
  if (hmm->flags & HMM_REF) strcpy(shmm->ref, hmm->ref);
  if (hmm->flags & HMM_CS)  strcpy(shmm->cs,  hmm->cs);
}

/* Function: PrintSearchHMM()
 * 
 * Purpose:  For debugging. Print out a search-style HMM
 */
void
PrintSearchHMM(FILE *fp, struct shmm_s *shmm)
{
  int k, x;
  int *tptr;

  fprintf(fp, " MATCH SCORES: \n");
  for (x = 0; x < 26; x++)
    {
      fprintf (fp, "%c ", 'A' + x);
      for (k = 0; k <= shmm->M; k++)
	fprintf(fp, "%5d ", shmm->m_emit[x][k]);
      fprintf(fp, "\n");
    }

  fprintf(fp, " INSERT SCORES: \n");
  for (x = 0; x < 26; x++)
    {
      fprintf (fp, "%c ", 'A' + x);
      for (k = 0; k <= shmm->M; k++)
	fprintf(fp, "%5d ", shmm->i_emit[x][k]);
      fprintf(fp, "\n");
    }

  printf("TRANSITION SCORES: \n");
  printf("     %5s %5s %5s %5s %5s %5s %5s %5s %5s \n",
	 "Tdd", "Tdi", "Tdm", 
	 "Tid", "Tii", "Tim", 
	 "Tmd", "Tmi", "Tmm");
  tptr = shmm->t;
  for (k = 0; k <= shmm->M; k++)
    {
      fprintf(fp, "     ");
      for (x = 0; x < 9; x++)
	{
	  fprintf(fp, "%5d ", *tptr);
	  tptr++;
	}
      fprintf(fp, "\n");
    }
}


struct shmm_s *
AllocSearchHMM(int M)
{
  struct shmm_s *shmm;
  int            x;

  if ((shmm = (struct shmm_s *) malloc (sizeof(struct shmm_s))) == NULL)
    Die("malloc failed");
  for (x = 0; x < 26; x++)
    if ((shmm->m_emit[x] = (int *) calloc (M+1, sizeof(int))) == NULL ||
	(shmm->i_emit[x] = (int *) calloc (M+1, sizeof(int))) == NULL)
      Die("malloc failed");
  if ((shmm->t   = (int *)  malloc (sizeof(int)  * (9*(M+1)))) == NULL ||
      (shmm->ref = (char *) malloc (sizeof(char) * (M+2))) == NULL ||
      (shmm->cs  = (char *) malloc (sizeof(char) * (M+2))) == NULL)
    Die("malloc failed");
  shmm->flags = 0;
  shmm->M = M;
  return shmm;
}

void
FreeSearchHMM(struct shmm_s *shmm)
{
  int x;

  for (x = 0; x < 26; x++)
    {
      free(shmm->m_emit[x]);
      free(shmm->i_emit[x]);
    }
  free(shmm->t);
  free(shmm->ref);
  free(shmm->cs);
  free(shmm);
}
  


/* Function: Symscore()
 * 
 * Purpose:  Given a sequence character x and an hmm containing
 *           probabilities, calculate the log-odds (base 2) score of
 *           the symbol for an emission scoring vector.
 *           
 * Args:     x      - the character, 'A'-'Z'
 *           scores - emission probability vector
 *           priors - prior probabilities for log-odds
 *                    
 * Return:   the integer log odds score of x given the emission
 *           vector and the priors, scaled up by INTSCALE.                   
 */
int
Symscore(char x, float *scores, float *priors)
{
  float  result;
  float  numer, denom;
  int    x_idx;

				/* simple case: x is in the alphabet */
  if (strchr(Alphabet, x) != NULL) 
    {
      x_idx  = SYMIDX(x);
      result = LOG2(scores[x_idx] / priors[x_idx]);
      return (int) (INTSCALE * result);
    }

  /* non-simple case: x is not in alphabet, but instead represents
   * an approved degenerate symbol (for instance, N for A|C|G|T.
   */
  if (Alphabet_type == kAmino)
    {
      switch (x) {
      case 'B': 
	numer  = scores[SYMIDX('N')] + scores[SYMIDX('D')];
	denom  = priors[SYMIDX('N')] + priors[SYMIDX('D')];
	break;
      case 'Z':
	numer  = scores[SYMIDX('Q')] + scores[SYMIDX('E')];
	denom  = priors[SYMIDX('Q')] + priors[SYMIDX('E')];
	break;
      default:
      case 'X':
	numer = denom = 1.0;
	break;
      }
    }
  else if (Alphabet_type == kDNA || Alphabet_type == kRNA)
    {
      switch (x) {		/* assumes order "ACGT" */
      case 'B':
	numer = scores[1] + scores[2] + scores[3]; 
	denom = priors[1] + priors[2] + priors[3]; 
	break;
      case 'D':	
	numer = scores[0] + scores[2] + scores[3];
	denom = priors[0] + priors[2] + priors[3];
	break;
      case 'H': 
	numer = scores[0] + scores[1] + scores[3];
	denom = priors[0] + priors[1] + priors[3];
	break;
      case 'K': 
	numer = scores[2] + scores[3];
	denom = priors[2] + priors[3];
	break;
      case 'M': 
	numer = scores[0] + scores[1];
	denom = priors[0] + priors[1];
	break;
      case 'R': 
	numer = scores[0] + scores[2];
	denom = priors[0] + priors[2];
	break;
      case 'S': 
	numer = scores[1] + scores[2];
	denom = priors[1] + priors[2];
	break;
      case 'V': 
	numer = scores[0] + scores[1] + scores[2];
	denom = priors[0] + priors[1] + priors[2];
	break;
      case 'T': 
      case 'U':
	numer = scores[3];
	denom = priors[3];
	break;
      case 'W':
	numer = scores[0] + scores[3];
	denom = priors[0] + priors[3];
	break;
      case 'Y':
	numer = scores[1] + scores[3];
	denom = priors[1] + priors[3];
	break;
      default:
      case 'N':
	numer = denom = 1.0;
	break;
      }
    }
  else
    {
      numer = denom = 1.0;
    }

  result = LOG2(numer / denom);
  return (INTSCALE * result);
}



/* Function: PrepareSequence()
 * 
 * Purpose:  Get a sequence ready for alignment to an HMM.
 *           Much like strdup, with the following exceptions:
 *              1) returns the length of the sequence, L
 *              2) string is converted to upper case
 *              3) any non-alphabet chars are converted to X
 *              4) the returned string is indexed 1..L, not 0..L-1
 *              5) there are dummy X's at 0 and L+1, for boundary
 *                 conditions; all the alignment algorithms
 *                 rely on these positions having zero emission
 *                 scores.
 *                 
 * Args:     s       - sequence to prepare
 *           ret_seq - RETURN: prepared sequence
 *           ret_L   - RETURN: length of sequence
 *           
 * Return:   (void)
 *           ret_seq is alloc'ed here and must be free'd by caller.          
 */                    
void
PrepareSequence(char *s, char **ret_seq, int *ret_L)
{
  int   L;
  int   pos;
  int   badchar = 0;
  char *seq;

  L = strlen(s);
  if ((seq = (char *) malloc (sizeof(char) * (L+3))) == NULL)
    Die("malloc failed");
  strcpy(seq+1, s);
  for (pos = 1; pos <= L; pos++)
    {
      if (islower(seq[pos])) 
	seq[pos] = toupper(seq[pos]);
      else if (! isalpha(seq[pos])) 
	{ seq[pos] = 'X'; badchar++; }
    }
  seq[0] = seq[L+1] = 'X';
  seq[L+2] = '\0';
  if (badchar) 
    fprintf(stderr, "Warning: converted %d non-alphabetic characters to X\n", badchar);
  *ret_seq = seq;
  if (ret_L != NULL) *ret_L   = L;
}


/* Function: ViterbiFill()
 *
 * Perform the Viterbi dynamic programming calculation step, in
 * its matrix fill stage. Columns are counted by k and range
 * 0..M+1; 1..M are the M positions of the HMM. Rows are counted
 * by i and range 0..L+1; 1..L are the L positions of the symbol
 * sequence.
 *
 * A non-allocated pointer "ret_mx" is passed; upon return, it will
 * point at the filled matrix. The caller is responsible for freeing
 * the allocated space.
 * 
 * The caller can obtain the score for the model aligned to the sequence
 * by getting the value in ret_mx[L+1][M+1].score_m.
 *
 * Return 0 on failure, 1 on success.
 */
int
ViterbiFill(struct shmm_s  *shmm,    /* integer log odds model       */
	    char           *seq,     /* sequence 1..L <- note offset */
	    int             L,       /* length of sequence           */
	    struct vit_s ***ret_mx,  /* RETURN: the calc'ed matrix   */
	    float          *ret_sc)  /* RETURN: best score           */
{
  struct vit_s **mx;            /* the viterbi calculation grid       */
  int    score;	                /* tmp variable for scores            */
  int    i;			/* counter for sequence position: 0,1..L */
  int    k;			/* counter for model position: 0,1..M */
  int    i_symidx;		/* hint for symbol index in alphabet  */
  struct vit_s *thisrow;        /* ptr to current row of mx           */
  struct vit_s *nextrow;        /* ptr to next row of mx              */
  int   *m_emit;                /* match scores for sym at this row   */
  int   *i_emit;                /* insert scores for sym at this row  */
  int   *tptr;                  /* transition scores                  */

  /********************************************
   * Initial setup and allocations
   ********************************************/
				/* allocate the calculation matrix,
				   which is 0..L+1 rows by 0..M+1 cols */
  if (( mx = (struct vit_s **) malloc (sizeof(struct vit_s *) * (L+2))) == NULL) 
      Die("memory failure allocating viterbi matrix\n");
  for (i = 0; i <= L+1; i++)
    if ((mx[i] = (struct vit_s *) malloc (sizeof(struct vit_s) * (shmm->M+2))) == NULL) 
      Die("memory failure allocating viterbi matrix, row %d\n", i);
  
  /********************************************
   * Initialization
   ********************************************/
      
				/* initialize the first cell 0,0 */
  mx[0][0].score_m = 0;
  mx[0][0].score_d = -99999999;
  mx[0][0].score_i = -99999999;

				/* initialize the top row */
  for (k = 1; k <= shmm->M+1; k++)
    {
      mx[0][k].score_m = -99999999;
      mx[0][k].score_i = -99999999;
    }

  /********************************************
   * Recursion: fill in the mx matrix
   ********************************************/

  for (i = 0; i <= L; i++)
    {
				/* get ptrs into current and next row. */
      thisrow = mx[i];
      nextrow = mx[i+1];
				/* initialize in the next row */
      nextrow[0].score_m = -99999999;
      nextrow[0].score_d = -99999999;


				/* setup scoring pointers for this row */
      i_symidx = seq[i] - 'A';
      m_emit = shmm->m_emit[i_symidx];
      i_emit = shmm->i_emit[i_symidx];
      tptr   = shmm->t;

      for (k = 0; k <= shmm->M; k++)
	{ /* begin inner loop... this is where all the time is spent. */

				/* add in emission scores to the current cell. */
	  if (i > 0)
	    {
	      thisrow[k].score_m += *m_emit; m_emit++;
	      thisrow[k].score_i += *i_emit; i_emit++;
	    }

				/* deal with transitions out of delete state */
				/* note: tptr assumes order dd,di,dm,id,ii,im,md,mi,mm */
				/* we initialize with these */
	  thisrow[k+1].score_d = thisrow[k].score_d + *tptr; tptr++;
	  nextrow[k].score_i   = thisrow[k].score_d + *tptr; tptr++;
	  nextrow[k+1].score_m = thisrow[k].score_d + *tptr; tptr++;

				/* deal with transitions out of insert state */
	  if ((score = thisrow[k].score_i + *tptr) > thisrow[k+1].score_d)
	    thisrow[k+1].score_d = score;
	  tptr++;
	  if ((score = thisrow[k].score_i + *tptr) > nextrow[k].score_i)
	    nextrow[k].score_i = score;
	  tptr++;
	  if ((score = thisrow[k].score_i + *tptr) > nextrow[k+1].score_m)
	    nextrow[k+1].score_m = score;
	  tptr++;

				/* deal with transitions out of match state */
	  if ((score = thisrow[k].score_m + *tptr) > thisrow[k+1].score_d)
	    thisrow[k+1].score_d = score;
	  tptr++;
	  if ((score = thisrow[k].score_m + *tptr) > nextrow[k].score_i)
	    nextrow[k].score_i = score;
	  tptr++;
	  if ((score = thisrow[k].score_m + *tptr) > nextrow[k+1].score_m)
	    nextrow[k+1].score_m = score;
	  tptr++;

	}
    }
  
  /********************************************
   * Termination of the last row (delete paths)
   ********************************************/
				/* only need to check transitions to deletes */
  tptr    = shmm->t;
  thisrow = mx[L+1];
  for (k = 0; k <= shmm->M; k++)
    {
      thisrow[k+1].score_d = thisrow[k].score_d + *tptr; 
      tptr += 3;
      if ((score = thisrow[k].score_i + *tptr) > thisrow[k+1].score_d)
	thisrow[k+1].score_d = score;
      tptr += 3;
      if ((score = thisrow[k].score_m + *tptr) > thisrow[k+1].score_d)
	thisrow[k+1].score_d = score;
      tptr += 3;
    }

  /********************************************
   * Return
   ********************************************/
  *ret_sc = (float) mx[L+1][shmm->M+1].score_m / INTSCALE;
  *ret_mx = mx;
  return 1;
}


/* Function: FreeViterbiMatrix()
 * 
 * Free a matrix filled by ViterbiFill().
 * 
 */
void
FreeViterbiMatrix(struct vit_s **mx,
		  int            L)  /* length of sequence matrix was for */
{
  int i;

  for (i = 0; i < L+2; i++)
    free(mx[i]);
  free(mx);
}


/* Function: ViterbiTrace()
 * 
 * Purpose:  General traceback function used by all but the simulated
 *           annealing variant of the alignment algorithms: VitFill(),
 *           SWViterbi, DBViterbi().
 *           
 *           The traceback rpos coordinates are 0..L-1.
 *           
 * Args:     mx     - filled viterbi scoring matrix
 *           shmm   - model, in scoring (integer log odds) form
 *           seq    - sequence, 1..L (note offset)
 *           window - number of rows (seq positions) that mx is valid for,
 *                    inclusive of from_i (Vit: L+2. DB,SW: window or end_i+1)
 *           end_i  - start tracing from this row (seq position, 0..L+1)
 *           end_k  - start tracing from this column (HMM pos, 0..M+1)
 *           ret_tr - RETURN: traceback structure 
 *           ret_i  - RETURN: starting row for alignment (seq position)
 *           ret_k  - RETURN: starting col for alignment (HMM position)
 *                    
 * Return: 1 on success, 0 if traceback exceeds window size and fails.
 *         ret_tr is alloc'ed here and must be freed by caller, using
 *         FreeTrace().
 */
int
ViterbiTrace(struct vit_s **mx, struct shmm_s *shmm, char *seq, int window,
	     int end_i, int end_k, 
	     struct trace_s **ret_tr, int *ret_i, int *ret_k)
{
  int  i;			/* counter for rows, 0..L+1 */
  int  curr;			/* counter for window position of row i */
  int  prev;
  int  k;			/* counter for cols, 0..M+1 */
  int  j;			/* counter for tmp state sequence, 0..L+M+1 */
  struct trace_s *tr;           /* traceback                                */
  int  N;			/* length of optimal state seq, for return */
  int  emitscore;

  /**************************************************
   * Allocation for traceback 
   **************************************************/
				/* we know N <= window+M+1, so alloc accordingly */
				/* leaving space for dummy BEGIN and END  */
  AllocTrace(window + shmm->M +3, &tr);

  /**************************************************
   * Initialization
   **************************************************/

  /* Dummy end state, and initial state, on traceback
   */
  j                = 0;
  tr->nodeidx[j]   = shmm->M+1;
  tr->statetype[j] = MATCH;
  tr->rpos[j]      = -1;
  i = end_i;       
  k = end_k;

  /* If end_k != shmm->M+1, we've got an internal match to
   * the model from SWViterbi(), and we need to add one
   * more initial state to the traceback.
   */
  if (end_k != shmm->M+1)
    {
      j++;
      tr->nodeidx[j]   = end_k;
      tr->statetype[j] = MATCH;
      tr->rpos[j]      = end_i-1;
    }

  /**************************************************
   * Traceback
   **************************************************/

  while (tr->statetype[j] != FROM_NOWHERE)
    {
      j++;
      prev = (i-1) % window;
      curr = i % window;

      if (end_i - i >= window)
	{			/* traceback exceeds window size. Quietly return failure. */
	  FreeTrace(tr);
	  return 0;
	}
				/* if we look back one position in tmp_statetype,
				   we find out what substate we're supposed to use
				   here in i,k */
      tr->statetype[j] = FROM_NOWHERE;
      switch (tr->statetype[j-1]) {
      case MATCH:		
	if (i > 0 && k > 0)
	  {
	    emitscore = (k == shmm->M+1)? 0 : shmm->m_emit[seq[i]-'A'][k];
	    
	    if (mx[prev][k-1].score_m + shmm->t[(k-1)*9+Tmm] + emitscore == mx[curr][k].score_m)
	      tr->statetype[j] = MATCH;
	    if (mx[prev][k-1].score_i + shmm->t[(k-1)*9+Tim] + emitscore == mx[curr][k].score_m)
	      tr->statetype[j] = INSERT;
	    if (mx[prev][k-1].score_d + shmm->t[(k-1)*9+Tdm] + emitscore == mx[curr][k].score_m)
	      tr->statetype[j] = DELETE;
	  }
				/* the state that got us here had index k-1, 
				   and we move to i-1, k-1 next */
	tr->nodeidx[j] = k-1;
	i--;
	k--;
	break;

      case DELETE:		
	if (k > 0)
	  {
	    if (mx[curr][k-1].score_m + shmm->t[(k-1)*9+Tmd] == mx[curr][k].score_d)
	      tr->statetype[j] = MATCH;
	    if (mx[curr][k-1].score_i + shmm->t[(k-1)*9+Tid] == mx[curr][k].score_d)
	      tr->statetype[j] = INSERT;
	    if (mx[curr][k-1].score_d + shmm->t[(k-1)*9+Tdd] == mx[curr][k].score_d)
	      tr->statetype[j] = DELETE;
	  }
				/* the state that got us here had index k-1, 
				   and we move to i, k-1 next */
	tr->nodeidx[j] = k-1;
	k--;
	break;

      case INSERT:		
	if (i > 0)
	  {
	    emitscore = (k == shmm->M+1)? 0 : shmm->i_emit[seq[i]-'A'][k];
	    
	    if (mx[prev][k].score_m + shmm->t[k*9+Tmi] + emitscore == mx[curr][k].score_i)
	      tr->statetype[j] = MATCH;
	    if (mx[prev][k].score_i + shmm->t[k*9+Tii] + emitscore == mx[curr][k].score_i)
	      tr->statetype[j] = INSERT;
	    if (mx[prev][k].score_d + shmm->t[k*9+Tdi] + emitscore == mx[curr][k].score_i)
	      tr->statetype[j] = DELETE;
	  }
				/* the state that got us here had index k, 
				   and we move to i-1, k next */
	tr->nodeidx[j] = k;
	i--;
	break;

      default:
	Die("unrecognized state type %d in traceback, i=%d k=%d j=%d\n",
	    tr->statetype[j-1], i, k, j);
      }

      if (tr->statetype[j] == INSERT || tr->statetype[j] == MATCH)
	tr->rpos[j] = i-1;
      else
	tr->rpos[j] = -1;
    }

  /* If we've run past the beginning of the model,
   * we did a complete match and have to just fix up
   * the last traceback pt we wrote; otherwise, we
   * must be doing Smith/Waterman, and we need to
   * add a completely extra cap on the trace
   */
  if (k < 0) j--;

  /* If we couldn't extend the trace any more, we're done.
   * Add a dummy begin state to cap it off.
   */
  tr->nodeidx[j]   = 0;
  tr->statetype[j] = MATCH;
  tr->rpos[j]      = -1;
  j++;

  /**************************************************
   * Reversing the state sequence
   **************************************************/
  
  N = j;			/* now we know the length of optimal state seq */
  ReverseTrace(tr, N);

  /**************************************************
   * Return
   **************************************************/

  /* Since we always end on a "match" state, we've always done a
   * i--, k-- as the last thing; so the true start position is i+1, k+1
   */
  if (ret_tr != NULL) *ret_tr = tr; else FreeTrace(tr);
  if (ret_i  != NULL) *ret_i  = i+1;
  if (ret_k  != NULL) *ret_k  = k+1;
  return 1;
}



/* Function: PrintViterbiMatrix()
 * 
 * Purpose:  Print out a (small) viterbi scoring matrix,
 *           for debugging purposes.
 *           
 * Args:     mx   - the score matrix, L+2 rows by M+2 columns
 *           seq1 - sequence (1..L)
 *           L    - length of seq
 *           M    - number of states in model                
 */
void
PrintViterbiMatrix(struct vit_s **mx, char *seq1, int L, int M)
{
  int i,k;

  for (i = 0; i <= L+1; i++)
    {
      printf("%d MAT: ",i);
      for (k = 0; k <= M+1; k++)
	printf("%10d ", mx[i][k].score_m);
      puts("");

      printf("%c INS: ", seq1[i]);
      for (k = 0; k <= M+1; k++)
	printf("%10d ", mx[i][k].score_i);
      puts("");
      
      printf("  DEL: ");
      for (k = 0; k <= M+1; k++)
	printf("%10d ", mx[i][k].score_d);
      puts("");
      puts("");
    }
}


/* Function: ViterbiAlignAlignment()
 * 
 * Purpose:  Align a multiple sequence alignment to an HMM without
 *           altering the multiple alignment.
 *           
 * Args:     shmm   - HMM in integer log-odds score form
 *           aseq   - alignment, [0..nseq-1][0..alen-1]
 *           alen   - length of aligned sequences
 *           nseq   - number of aligned sequences
 *           ret_tr - RETURN: array of tracebacks. rpos field is
 *                    relative to aseq, not raw seq, similar to
 *                    Maxmodelmaker(); use DealignTrace() if you
 *                    want relative to raw sequence.
 *           ret_sc - RETURN: sum of log odds scores.
 *           
 * Return:   (void)
 *           ret_tr is alloced here. Individuals must be free'd by FreeTrace(),
 *           then tr itself free'd by free().
 */
void
ViterbiAlignAlignment(struct shmm_s *shmm, char **aseq, int alen, int nseq,
		      struct trace_s ***ret_tr, float *ret_sc)
{
  struct fvit_s **mx;           /* the viterbi calculation grid       */
  int    score;	                /* tmp variable for scores            */
  int    i;			/* counter for sequence position: 0,1..L */
  int    k;			/* counter for model position: 0,1..M */
  int    idx;			/* index for sequences                */
  struct fvit_s *thisrow;       /* ptr to current row of mx           */
  struct fvit_s *nextrow;       /* ptr to next row of mx              */
  int  **matocc;                /* [0..alen+1][0..nseq-1], 1 for MATCH*/
  struct trace_s **tr;          /* array of tracebacks to return      */
  int   *tpos;                  /* index for position in indiv traces */
  int    lastsub;		/* last state type in master trace    */

  /* A crucial extra component of this alignment algorithm:
   * at each matrix cell, we have to remember: for the best
   * path into the INSERT subcell, what state is each sequence in?
   * This is non-trivial because some gaps are assigned to
   * no states. When we calculate the score from an insert column,
   * where there are gaps we have to look up the previous state.
   *
   * Fortunately, we don't need to keep a full matrix of these,
   * or we'd be in serious memory problems. Use a rolling pointer
   * trick, keep two active rows "current" and "next".
   */
  char **cur_state;             /* [0..M+1][0..nseq-1]; MATCH, INSERT, or DELETE */ 
  char **nxt_state;             /* same, except keeps states for next row        */
  char **swap;                  /* used for swapping cur, nxt                    */

  /********************************************
   * Initial setup and allocations
   ********************************************/
				/* allocate the calculation matrix,
				   which is 0..alen+1 rows by 0..M+1 cols */
  mx        = (struct fvit_s **) MallocOrDie (sizeof(struct fvit_s *) * (alen+2));
  matocc    = (int  **)          MallocOrDie (sizeof(int *)           * (alen+2));
  cur_state = (char **)          MallocOrDie (sizeof(char *)          * (shmm->M+2));
  nxt_state = (char **)          MallocOrDie (sizeof(char *)          * (shmm->M+2));
  for (i = 0; i <= alen+1; i++)
    {
      mx[i]    = (struct fvit_s *) MallocOrDie (sizeof(struct fvit_s) * (shmm->M+2));
      matocc[i]= (int *)           MallocOrDie (sizeof(int)           * nseq);
    }
  for (k = 0; k <= shmm->M+1; k++)
    {
      cur_state[k] = (char *) MallocOrDie (sizeof(char) * nseq);
      nxt_state[k] = (char *) MallocOrDie (sizeof(char) * nseq);
    }

  /********************************************
   * Initialization
   ********************************************/
				/* initialize the first cell 0,0 */
  mx[0][0].score_m = 0;
  mx[0][0].score_d = -99999999;
  mx[0][0].score_i = -99999999;

  for (k = 0; k <= shmm->M+1; k++)
    for (idx = 0; idx < nseq; idx++)
      nxt_state[k][idx] = MATCH;

				/* initialize the top row */
  for (k = 1; k <= shmm->M+1; k++)
    {
      mx[0][k].score_m = -99999999;
      mx[0][k].score_i = -99999999;
    }

  /* Precalculate matocc (match occupancy). 
   * 1 if symbol in column for this seq, 0 if not. 
   * 1..alen, from 0..alen-1 alignments
   */
  for (idx = 0; idx < nseq; idx++)
    {
      matocc[0][idx] = matocc[alen+1][idx] = 1; /* dummies for BEGIN, END */
      for (i = 1; i <= alen; i++)
	matocc[i][idx] = isgap(aseq[idx][i-1]) ? 0 : 1;
    }

  /********************************************
   * Recursion: fill in the mx matrix
   ********************************************/
				/* Alignment is 0..alen-1, we index it 
				   here as 1..alen because of Viterbi matrix. */
  for (i = 0; i <= alen; i++)
    {
				/* get ptrs into current and next row. */
      thisrow = mx[i];
      nextrow = mx[i+1];
				/* initialize in the next row */
      nextrow[0].score_m = -99999999;
      nextrow[0].score_d = -99999999;

      swap = cur_state; cur_state = nxt_state; nxt_state = swap;

      for (k = 0; k <= shmm->M; k++)
	{ /* begin inner loop... this is where all the time is spent. */

				/* add in emission scores to the current cell. */
	  if (i > 0)
	    for (idx = 0; idx < nseq; idx++)
	      if (matocc[i][idx])
		{
		  thisrow[k].score_m += shmm->m_emit[aseq[idx][i-1] - 'A'][k];
		  thisrow[k].score_i += shmm->i_emit[aseq[idx][i-1] - 'A'][k];
		}
				/* initialize with transitions out of delete state */
				/* to delete */
	  thisrow[k+1].score_d = thisrow[k].score_d + shmm->t[9*k + Tdd] * nseq;
	  thisrow[k+1].tback_d = DELETE;
				/* to insert */
	  nextrow[k].score_i = thisrow[k].score_d;
	  nextrow[k].tback_i = DELETE;
	  for (idx = 0; idx < nseq; idx++) 
	    if (matocc[i+1][idx])
	      {
		nextrow[k].score_i += shmm->t[9*k + Tdi];
		nxt_state[k][idx]  = INSERT;
	      }
	    else
	      nxt_state[k][idx] = DELETE;
				/* to match */
	  nextrow[k+1].score_m = thisrow[k].score_d;
	  nextrow[k+1].tback_m = DELETE;
	  for (idx = 0; idx < nseq; idx++)
	    if (matocc[i+1][idx])
	      nextrow[k+1].score_m += shmm-> t[9*k + Tdm];
	    else
	      nextrow[k+1].score_m += shmm-> t[9*k + Tdd];

	  
				/* deal with transitions out of insert state */
				/* to delete state. */
	  score = thisrow[k].score_i;
	  for (idx = 0; idx < nseq; idx++)
	    switch (cur_state[k][idx]) {
	    case MATCH:  score += shmm->t[9*k + Tmd]; break;
	    case DELETE: score += shmm->t[9*k + Tdd]; break;
	    case INSERT: score += shmm->t[9*k + Tid]; break;
	    }
	  if (score > thisrow[k+1].score_d) 
	    {
	      thisrow[k+1].score_d = score;
	      thisrow[k+1].tback_d = INSERT;
	    }
				/* to insert state */
	  score = thisrow[k].score_i;
	  for (idx = 0; idx < nseq; idx++)
	    {
	      if (matocc[i+1][idx])
		switch (cur_state[k][idx]) {
		case MATCH:  score += shmm->t[9*k + Tmi]; break;
		case DELETE: score += shmm->t[9*k + Tdi]; break;
		case INSERT: score += shmm->t[9*k + Tii]; break;
		}
	    }
	  if (score > nextrow[k].score_i) 
	    {
	      nextrow[k].score_i = score;
	      nextrow[k].tback_i = INSERT;
	      for (idx = 0; idx < nseq; idx++)
		if (matocc[i+1][idx])
		  nxt_state[k][idx] = INSERT;
		else
		  nxt_state[k][idx] = cur_state[k][idx];
	    }
				/* to match state */
	  score = thisrow[k].score_i;
	  for (idx = 0; idx < nseq; idx++)
	    if (matocc[i+1][idx])
	      switch (cur_state[k][idx]) {
	      case MATCH:  score += shmm->t[9*k + Tmm]; break;
	      case DELETE: score += shmm->t[9*k + Tdm]; break;
	      case INSERT: score += shmm->t[9*k + Tim]; break;
	      }
	    else
	      switch (cur_state[k][idx]) {
	      case MATCH:  score += shmm->t[9*k + Tmd]; break;
	      case DELETE: score += shmm->t[9*k + Tdd]; break;
	      case INSERT: score += shmm->t[9*k + Tid]; break;
	      }
	  if (score > nextrow[k+1].score_m) 
	    {
	      nextrow[k+1].score_m = score;
	      nextrow[k+1].tback_m = INSERT;
	    }

	  /* Transitions out of match state.
	   */
				/* to delete */
	  score = thisrow[k].score_m;
	  for (idx = 0; idx < nseq; idx++)
	    if (matocc[i][idx])
	      score += shmm->t[9*k + Tmd];
	    else
	      score += shmm->t[9*k + Tdd];
	  if (score > thisrow[k+1].score_d)
	    {
	      thisrow[k+1].score_d = score;
	      thisrow[k+1].tback_d = MATCH;
	    }
				/* to insert */
	  score = thisrow[k].score_m;
	  for (idx = 0; idx < nseq; idx++)
	    if (matocc[i+1][idx])
	      {
		if (matocc[i][idx])
		  score += shmm->t[9*k + Tmi];
		else
		  score += shmm->t[9*k + Tdi];
	      }
	  if (score > nextrow[k].score_i)
	    {
	      nextrow[k].score_i = score;
	      nextrow[k].tback_i = MATCH;
	      for (idx = 0; idx < nseq; idx++)
		if (matocc[i+1][idx])
		  nxt_state[k][idx] = INSERT;
		else if (matocc[i][idx])
		  nxt_state[k][idx] = MATCH;
		else
		  nxt_state[k][idx] = DELETE;
	    }
				/* to match */
	  score = thisrow[k].score_m;
	  for (idx = 0; idx < nseq; idx++)
	    if (matocc[i][idx])
	      {
		if (matocc[i+1][idx])
		  score += shmm->t[9*k + Tmm];
		else
		  score += shmm->t[9*k + Tmd];
	      }
	    else
	      {
		if (matocc[i+1][idx])
		  score += shmm->t[9*k + Tdm];
		else
		  score += shmm->t[9*k + Tdd];
	      }
	  if (score > nextrow[k+1].score_m)
	    {
	      nextrow[k+1].score_m = score;
	      nextrow[k+1].tback_m = MATCH;
	    }

	} /* end loop over model positions k */

    } /* end loop over alignment positions i */

/*  PrintFragViterbiMatrix(mx, alen, shmm->M); */

  /* Fill stage finished.
   * mx now contains final score in mx[alen+1][M+1].
   * Trace back from there to get master alignment.
   */
  tr   = (struct trace_s **) MallocOrDie (sizeof(struct trace_s *) * nseq);
  tpos = (int *)             MallocOrDie (sizeof(int)              * nseq);
  for (idx = 0; idx < nseq; idx++)
    {
      AllocTrace(alen + shmm->M + 3, &(tr[idx]));
      tr[idx]->nodeidx[0]   = shmm->M+1;
      tr[idx]->statetype[0] = MATCH;
      tr[idx]->rpos[0]      = -1;
      tpos[idx]        = 1;
    }
  i      = alen+1;
  k      = shmm->M+1;
  lastsub= MATCH;
	       
  while (i != 0 || k != 0)
    {
      switch (lastsub) {
      case MATCH:  lastsub = mx[i][k].tback_m; i--; k--; break;
      case DELETE: lastsub = mx[i][k].tback_d;      k--; break;
      case INSERT: lastsub = mx[i][k].tback_i; i--;      break;
      default: Die("trace failed!");
      }

      switch (lastsub) {
      case MATCH:
	for (idx = 0; idx < nseq; idx++)
	  if (matocc[i][idx]) 
	    {
	      tr[idx]->nodeidx[tpos[idx]]   = k;
	      tr[idx]->statetype[tpos[idx]] = MATCH;
	      tr[idx]->rpos[tpos[idx]]      = i-1;
	      tpos[idx]++;
	    }
	  else
	    {
	      tr[idx]->nodeidx[tpos[idx]]   = k;
	      tr[idx]->statetype[tpos[idx]] = DELETE;
	      tr[idx]->rpos[tpos[idx]]      = -1;
	      tpos[idx]++;
	    }
	break;
      case INSERT:
	for (idx = 0; idx < nseq; idx++)
	  if (matocc[i][idx])
	    {
	      tr[idx]->nodeidx[tpos[idx]]   = k;
	      tr[idx]->statetype[tpos[idx]] = INSERT;
	      tr[idx]->rpos[tpos[idx]]      = i-1;
	      tpos[idx]++;
	    }
	break;
      case DELETE:
	for (idx = 0; idx < nseq; idx++)
	  {
	    tr[idx]->nodeidx[tpos[idx]]   = k;
	    tr[idx]->statetype[tpos[idx]] = DELETE;
	    tr[idx]->rpos[tpos[idx]]      = -1;
	    tpos[idx]++;
	  }
	break;
      default: Die("trace failed!");
      }	/* end switch across new subcell in traceback */
    } /* end traceback */

  for (idx = 0; idx < nseq; idx++)
    ReverseTrace(tr[idx], tpos[idx]);

  *ret_tr = tr;
  *ret_sc = (float) mx[alen+1][shmm->M+1].score_m / INTSCALE;

  Free2DArray(matocc, alen+2);
  Free2DArray(cur_state, shmm->M+2);
  Free2DArray(nxt_state, shmm->M+2);
  Free2DArray(mx, alen+2);
  free(tpos);
}


/* Function: PrintFragViterbiMatrix()
 * 
 * Purpose:  Print out a (small) viterbi scoring matrix,
 *           for debugging purposes.
 *           
 * Args:     mx   - the score matrix, L+2 rows by M+2 columns
 *           L    - length of seq
 *           M    - number of states in model                
 */
void
PrintFragViterbiMatrix(struct fvit_s **mx, int L, int M)
{
  int i,k;

  for (i = 0; i <= L+1; i++)
    {
      printf("%d MAT: ",i);
      for (k = 0; k <= M+1; k++)
	printf("%10d ", mx[i][k].score_m);
      puts("");

      printf("  INS: ");
      for (k = 0; k <= M+1; k++)
	printf("%10d ", mx[i][k].score_i);
      puts("");
      
      printf("  DEL: ");
      for (k = 0; k <= M+1; k++)
	printf("%10d ", mx[i][k].score_d);
      puts("");
      puts("");
    }
}

