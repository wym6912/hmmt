/* SQUID - A C function library for biological sequence analysis
 * Copyright (C) 1992-1995 Sean R. Eddy	
 *
 *    This source code is distributed under terms of the
 *    GNU General Public License. See the files COPYING 
 *    and GNULICENSE for further details.
 *
 */

/* weight.c
 * SRE, Thu Mar  3 07:56:01 1994
 * 
 * Calculate weights for sequences in an alignment.
 */

#include <ctype.h>
#include <string.h>
#include "squid.h"

static void upweight(struct phylo_s *tree, int nseq, float *lwt, float *rwt, int node);
static void downweight(struct phylo_s *tree, int nseq, float *lwt, float *rwt, 
		       float *fwt, int node);
static float simple_distance(char *s1, char *s2);
static int    simple_diffmx(char **aseqs,int num, float ***ret_dmx);

/* Function: SonnhammerWeights()
 * 
 * Purpose:  Use Erik's tree-based algorithm to set weights for
 *           sequences in an alignment. upweight() and downweight()
 *           are derived from Graeme Mitchison's code.
 *           
 * Args:     aseq        - array of (0..nseq-1) aligned sequences, all of 
 *                           same length
 *           nseq        - number of sequences in aseq
 *           alen        - lengths of all sequences (all must be the same length)
 *           ret_weights - RETURN: 0..nseq-1 sequence weights. Weights
 *                           sum to nseq.          
 *                           
 * Return:   1 on success, 0 on failure.
 *           caller is responsible for freeing *ret_weights 
 */
int
SonnhammerWeights(char **aseq, int nseq, int alen, float **ret_weights)
{
  float **dmx;                 /* distance (difference) matrix */
  struct phylo_s *tree;
  float  *lwt, *rwt;           /* weight on left, right of this tree node */
  float  *fwt;                 /* final weight assigned to this node */
  float  *wt;                  /* RETURN: sequence weights           */
  int      i;
  
  /* I use a simple fractional difference matrix derived by
   * pairwise identity. Perhaps I should include a Poisson
   * distance correction.
   */
  if (! MakeDiffMx(aseq, nseq, alen, TRUE, &dmx)) Die("MakeDiffMx() failed");
  if (! Cluster(dmx, nseq, CLUSTER_MIN, &tree))  Die("Cluster() failed");
  
  /* Allocations
   */
  if ((lwt = (float *) malloc (sizeof(float) * (2 * nseq - 1))) == NULL ||
      (rwt = (float *) malloc (sizeof(float) * (2 * nseq - 1))) == NULL ||
      (fwt = (float *) malloc (sizeof(float) * (2 * nseq - 1))) == NULL ||
      (wt  = (float *) malloc (sizeof(float) * (nseq)))         == NULL)
    Die("malloc failed");
  
  /* lwt and rwt are the total branch weight to the left and
   * right of a node or sequence. They are 0..2N-2.  0..N-1 are 
   * the sequences; these have weight 0. N..2N-2 are the actual
   * tree nodes.
   */
  for (i = 0; i < nseq; i++)
    lwt[i] = rwt[i] = 0.0;
				/* recursively calculate rwt, lwt, starting
				   at node nseq (the root) */
  upweight(tree, nseq, lwt, rwt, nseq);

				/* recursively distribute weight across the
				   tree */
  fwt[nseq] = nseq;
  downweight(tree, nseq, lwt, rwt, fwt, nseq);
				/* collect the weights */
  for (i = 0; i < nseq; i++)
    wt[i] = fwt[i];

  Free2DArray(dmx, nseq);
  FreePhylo(tree, nseq);
  free(lwt); free(rwt); free(fwt);
  *ret_weights = wt;
  return 1;
}

static void 
upweight(struct phylo_s *tree, int nseq, float *lwt, float *rwt, int node)
{
  int ld,rd;

  ld = tree[node-nseq].left;
  if (ld >= nseq) upweight(tree, nseq, lwt, rwt, ld);
  rd = tree[node-nseq].right;
  if (rd >= nseq) upweight(tree, nseq, lwt, rwt, rd);
  lwt[node] = lwt[ld] + rwt[ld] + tree[node-nseq].lblen;
  rwt[node] = lwt[rd] + rwt[rd] + tree[node-nseq].rblen;
}


static void 
downweight(struct phylo_s *tree, int nseq, float *lwt, float *rwt, float *fwt, int node)
{
  int ld,rd;
  float lnum, rnum;

  ld = tree[node-nseq].left;
  rd = tree[node-nseq].right;
  if (lwt[node] + rwt[node] > 0.0)
    {
      fwt[ld] = fwt[node] * (lwt[node] / (lwt[node] + rwt[node]));
      fwt[rd] = fwt[node] * (rwt[node] / (lwt[node] + rwt[node]));
    }
  else
    {
      lnum = (ld >= nseq) ? tree[ld-nseq].incnum : 1.0;
      rnum = (rd >= nseq) ? tree[rd-nseq].incnum : 1.0;
      fwt[ld] = fwt[node] * lnum / (lnum + rnum);
      fwt[rd] = fwt[node] * rnum / (lnum + rnum);
    }

  if (ld >= nseq) downweight(tree, nseq, lwt, rwt, fwt, ld);
  if (rd >= nseq) downweight(tree, nseq, lwt, rwt, fwt, rd);
}




/* Function: VoronoiWeights()
 * 
 * Purpose:  Calculate weights using the scheme of Sibbald &
 *           Argos (JMB 216:813-818 1990). The scheme is
 *           slightly modified because the original algorithm
 *           actually doesn't work on gapped alignments.
 *           The sequences are assumed to be protein.
 *           
 * Args:     aseq        - array of (0..nseq-1) aligned sequences, all of 
 *                           same length
 *           nseq        - number of sequences in aseq
 *           alen        - lengths of sequences (all the same length)
 *           ret_weights - RETURN: 0..nseq-1 sequence weights. Weights
 *                           sum to nseq.          
 *
 * Return:   1 on success, 0 on failure.
 *           caller is responsible for freeing *ret_weights 
 */
int
VoronoiWeights(char **aseq, int nseq, int alen, float **ret_weights)
{
  float **dmx;                 /* distance (difference) matrix    */
  float  *halfmin;             /* 1/2 minimum distance to other seqs */
  float  *wt;                  /* RETURN: sequence weights        */
  char   **psym;                /* symbols seen in each column     */
  int     *nsym;                /* # syms seen in each column      */
  int      symseen[27];         /* flags for observed syms         */
  char    *randseq;             /* randomly generated sequence     */
  int      acol;		/* pos in aligned columns          */
  int      idx;                 /* index in sequences              */
  int      symidx;              /* 0..25 index for symbol          */
  int      i;			/* generic counter                 */
  float   min;			/* minimum distance                */
  float   dist;		/* distance between random and real */
  long     challenge, champion; /* for resolving ties              */
  int      itscale;		/* how many iterations per seq     */
  int      iteration;           
  int      best;		/* index of nearest real sequence  */

  itscale = 50;

  /* Precalculate 1/2 minimum distance to other
   * sequences for each sequence
   */
  if (! simple_diffmx(aseq, nseq, &dmx)) 
    Die("simple_diffmx() failed");
  if ((halfmin = (float *) malloc (sizeof(float) * nseq)) == NULL)
    Die("malloc failed");
  for (idx = 0; idx < nseq; idx++)
    {
      for (min = 1.0, i = 0; i < nseq; i++)
	{
	  if (i == idx) continue;
	  if (dmx[idx][i] < min) min = dmx[idx][i];
	}
      halfmin[idx] = min / 2.0;
    }
  Free2DArray(dmx, nseq);

  /* Set up the random sequence generating model.
   */
  if ((psym = (char **) malloc (alen * sizeof(char *))) == NULL ||
      (nsym = (int *)   malloc (alen * sizeof(int)))    == NULL)
    Die("malloc failed");
  for (acol = 0; acol < alen; acol++)
    if ((psym[acol] = (char *) malloc (27 * sizeof(char))) == NULL)
      Die("malloc failed");

/* #ifdef ORIGINAL_SIBBALD_ALGORITHM_IS_BROKEN */
  for (acol = 0; acol < alen; acol++)
    {
      memset(symseen, 0, sizeof(int) * 27);
      for (idx = 0; idx < nseq; idx++)
	if (! isgap(aseq[idx][acol]))
	  {
	    if (isupper(aseq[idx][acol])) 
	      symidx = aseq[idx][acol] - 'A';
	    else
	      symidx = aseq[idx][acol] - 'a';
	    if (symidx >= 0 && symidx < 26)
	      symseen[symidx] = 1;
	  }
	else
	  symseen[26] = 1;	/* a gap */

      for (nsym[acol] = 0, i = 0; i < 26; i++)
	if (symseen[i]) 
	  {
	    psym[acol][nsym[acol]] = 'A'+i;
	    nsym[acol]++;
	  }
      if (symseen[26]) { psym[acol][nsym[acol]] = ' '; nsym[acol]++; }
    }
/* #endif ORIGINAL_SIBBALD_ALGORITHM_IS_BROKEN */

  /* Note: the original Sibbald&Argos algorithm calls for
   * bounding the sampled space using a template-like random
   * sequence generator. However, this leads to one minor
   * and one major problem. The minor problem is that
   * exceptional amino acids in a column can have a
   * significant effect by altering the amount of sampled
   * sequence space; the larger the data set, the worse
   * this problem becomes. The major problem is that 
   * there is no reasonable way to deal with gaps.
   * Gapped sequences simply inhabit a different dimensionality
   * and it's pretty painful to imagine calculating Voronoi
   * volumes when the N in your N-space is varying.
   * Note that all the examples shown by Sibbald and Argos
   * are *ungapped* examples.
   * 
   * The best way I've found to circumvent this problem is
   * just not to bound the sampled space; count gaps as
   * symbols and generate completely random sequences.
   */
#ifdef ALL_SEQUENCE_SPACE
  for (acol = 0; acol < alen; acol++)
    {
      strcpy(psym[acol], "ACDEFGHIKLMNPQRSTVWY ");
      nsym[acol] = 21;
    }
#endif
  
  /* Sibbald and Argos algorithm:
   *   1) assign all seqs weight 0.
   *   2) generate a "random" sequence
   *   3) calculate distance to every other sequence
   *      (if we get a distance < 1/2 minimum distance
   *       to other real seqs, we can stop)
   *   4) if unique closest sequence, increment its weight 1.
   *      if multiple closest seq, choose one randomly    
   *   5) repeat 2-4 for lots of iterations
   *   6) normalize all weights to sum to nseq.
   */
  if ((randseq = (char *) malloc ((alen+1) * sizeof(char))) == NULL ||
      (wt = (float *) malloc (nseq * sizeof(float))) == NULL)
    Die("malloc failed");
  for (idx = 0; idx < nseq; idx++)
    wt[idx] = 0.0;
  for (iteration = 0; iteration < itscale * nseq; iteration++)
    {
      for (acol = 0; acol < alen; acol++)
	randseq[acol] = (nsym[acol] == 0) ? ' ' : psym[acol][CHOOSE(nsym[acol])];
      randseq[acol] = '\0';

      champion = sre_random();
      for (min = 1.0, idx = 0; idx < nseq; idx++)
	{
	  dist = simple_distance(aseq[idx], randseq);
	  if (dist < halfmin[idx]) 
	    { 
	      best = idx; 
	      break;      
	    } 
	  if (dist < min)          
	    { champion = sre_random(); best = idx; min = dist; } 
	  else if (dist == min)    
	    { 
	      challenge = sre_random(); 
	      if (challenge > champion)
		{ champion = challenge; best = idx; min = dist; }
	    }
	}
      wt[best] += 1.0;
    }
  for (idx = 0; idx < nseq; idx++)
    wt[idx] = wt[idx] / (float) itscale;

  free(randseq);
  free(nsym);
  free(halfmin);
  Free2DArray(psym, alen);
  *ret_weights = wt;
  return 1;
}


/* Function: simple_distance()
 * 
 * Purpose:  For two identical-length null-terminated strings, return
 *           the fractional difference between them. (0..1)
 *           (Gaps don't count toward anything.)
 */
static float
simple_distance(char *s1, char *s2)
{
  int diff  = 0;
  int valid = 0;

  for (; *s1 != '\0'; s1++, s2++)
    {
      if (isgap(*s1) || isgap(*s2)) continue;
      if (*s1 != *s2) diff++;
      valid++;
    }
  return (float) diff / (float) valid;
}
    
/* Function: simple_diffmx()
 * 
 * Purpose:  Given a set of flushed, aligned sequences, construct
 *           an NxN fractional difference matrix using the
 *           simple_distance rule.
 *           
 * Args:     aseqs        - flushed, aligned sequences
 *           num          - number of aseqs
 *           ret_dmx      - RETURN: difference matrix (caller must free)
 *           
 * Return:   1 on success, 0 on failure.
 */
static int
simple_diffmx(char    **aseqs,
	      int       num,
	      float ***ret_dmx)
{
  float **dmx;                 /* RETURN: distance matrix           */
  int      i,j;			/* counters over sequences           */

  /* Allocate
   */
  if ((dmx = (float **) malloc (sizeof(float *) * num)) == NULL)
    Die("malloc failed");
  for (i = 0; i < num; i++)
    if ((dmx[i] = (float *) malloc (sizeof(float) * num)) == NULL)
      Die("malloc failed");

  /* Calculate distances, symmetric matrix
   */
  for (i = 0; i < num; i++)
    for (j = i; j < num; j++)
      dmx[i][j] = dmx[j][i] = simple_distance(aseqs[i], aseqs[j]);

  /* Return
   */
  *ret_dmx = dmx;
  return 1;
}
