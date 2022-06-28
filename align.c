/************************************************************
 * HMMER - Biological sequence analysis with HMMs
 * Copyright 1992-1995 Sean R. Eddy
 *
 *   This source code is distributed under the terms of the
 *   GNU General Public License. See the files COPYING and
 *   GNULICENSE for details.
 *
 ************************************************************/

/* align.c
 * Multiple sequence alignment generation, from raw sequences
 * and their tracebacks.
 * 
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#include "states.h"
#include "externs.h"
#include "squid.h"
#include "version.h"

#ifdef MEMDEBUG
#include "dbmalloc.h"
#endif


/* Function: Traces2Alignment()
 * 
 * Purpose:  Given tracebacks of some sequences rseq, construct
 *           and return a (flushed) multiple sequence alignment.
 *           
 * Args:     rseqs      - raw sequences being aligned          
 *           sqinfo     - array of info about the sequences
 *           nseq       - number of rseqs
 *           M          - length of model
 *           tr         - array of tracebacks
 *           matchonly  - TRUE if we don't print insert-generated symbols at all
 *           ret_aseqs  - RETURN: multiple sequence alignment           
 *           ainfo      - RETURN: info about alignment (alen, structure, and sqinfo)
 *           
 * Return:   1 on success, 0 on failure.
 */          
int
Traces2Alignment(char **rseqs, SQINFO *sqinfo, int nseq, int M, 
		 struct trace_s **tr, int matchonly, 
		 char ***ret_aseqs, AINFO *ainfo)
{
  char **aseqs;                 /* RETURN: aligned sequence set */
  int    seqidx;		/* counter for sequences */
  int    alen;		        /* width of alignment */
  int   *inserts;               /* array of maximum gaps between aligned columns */
  int    count;			/* counter for inserts */
  int    apos;			/* position in aligned sequence */
  int    rpos;			/* position in raw sequence */
  int    statepos;		/* position counter in state sequence (path) */
  int    k;			/* counter over states in model */

  /* Here's the problem. We want to align the match states in columns,
   * but some sequences have inserted symbols in them; we need some
   * sort of overall knowledge of where the inserts are and how long
   * they are in order to create the alignment.
   * 
   * Here's our trick. inserts[] is a 0..hmm->M array; inserts[i] stores
   * the maximum number of times insert substate i was used. This
   * is the maximum number of gaps to insert between canonical 
   * column i and i+1. 
   */

  if ((inserts = (int *) malloc (sizeof(int) * (M+1))) == NULL)
    Die("malloc failed for inserts array");

  for (k = 0; k <= M; k++)
    inserts[k] = 0;

  for (seqidx = 0; seqidx < nseq; seqidx++)
    {
      count = 0;
      for (statepos = 1; statepos < tr[seqidx]->tlen; statepos++)
	{
	  switch (tr[seqidx]->statetype[statepos]) 
	    {
	    case MATCH:
	    case DELETE:
	      if (inserts[tr[seqidx]->nodeidx[statepos]-1] < count)
		inserts[tr[seqidx]->nodeidx[statepos]-1] = count;
	      count = 0;
	      break;
	    case INSERT:
	      count++;
	      break;
	    default:
	      Die("Traces2Alignment reports unrecognized statetype %c",
		  tr[seqidx]->statetype[statepos]);
	    }
	}
    }

  /***********************************************
   * Construct the alignment
   ***********************************************/
				/* figure out the largest possible aligned string */
  alen = M;
  for (k = 0; k <= M ; k++)
    if (!matchonly) alen += inserts[k];
    else if (inserts[k] > 0) alen++;

  				/* allocations */
  if ((aseqs = (char **) malloc (sizeof(char *) * nseq)) == NULL)
    Die("malloc failed for aseqs array");

  for (seqidx = 0; seqidx < nseq; seqidx++)
    if ((aseqs[seqidx] = (char *) malloc (alen + 1)) == NULL)
      Die("malloc failed for aseqs member %d\n", seqidx);
  
				/* write each sequence */
  for (seqidx = 0; seqidx < nseq; seqidx ++)
    {
      apos = 0;
      rpos = 0;

				/* watch out for fragment alignments */
      if (tr[seqidx]->nodeidx[1] > 1)
	for (k = 1; k < tr[seqidx]->nodeidx[1]; k++)
	  {
	    if (! matchonly)	/* deal with insert states */
	      for (count = 0; count < inserts[k-1]; count++)
		aseqs[seqidx][apos++] = '.';
	    aseqs[seqidx][apos++] = '.'; /* the match state skipped by fragment */
	  }
				/* avoid BEGIN and END state */
      count = 0;
      for (statepos = 1; statepos < tr[seqidx]->tlen-1; statepos++)
	{
	  switch (tr[seqidx]->statetype[statepos]) 
	    {
	    case MATCH:
				/* make up gaps */
	      if (! matchonly)
		{		/* Right-justify half of the insertion  */
		  int rightshove, iold, inew;
		  char sym;

		  rightshove = count / 2;
		  iold       = apos - 1;
		  for (; count < inserts[tr[seqidx]->nodeidx[statepos]-1]; count++)
		    {
		      aseqs[seqidx][apos] = '.';
		      apos++;
		    }
		  inew = apos-1;
		  for (count = 0; count < rightshove; count++, iold--, inew--)
		    {
		      sym = aseqs[seqidx][iold];
		      aseqs[seqidx][iold] = '.';
		      aseqs[seqidx][inew] = sym;
		    }
		}
				/* if we're only showing match states,
				   we put a single space where at least one
				   insertion showed up */
	      else if (inserts[tr[seqidx]->nodeidx[statepos]-1] > 0)
		{
		  aseqs[seqidx][apos] = count ? 'x' : '.';
		  apos++;
		}
				/* emit the match */
	      aseqs[seqidx][apos] = rseqs[seqidx][rpos];
	      apos++;
	      rpos++;
	      count = 0;
	      break;

	    case DELETE:
				/* make up gaps */
	      if (! matchonly)
		for (; count < inserts[tr[seqidx]->nodeidx[statepos]-1]; count++)
		  {
		    aseqs[seqidx][apos] = '.';
		    apos++;
		  }
	      else if (inserts[tr[seqidx]->nodeidx[statepos]-1] > 0)
		{
		  aseqs[seqidx][apos] = count ? 'x' : '.';
		  apos++;
		}
		
				/* emit the deletion */
	      aseqs[seqidx][apos] = '.';
	      apos++;
	      count = 0;
	      break;
	    case INSERT:
	      if (! matchonly)
		{
		  aseqs[seqidx][apos] = sre_tolower((int) rseqs[seqidx][rpos]);
		  apos++;
		}
	      count++;
	      rpos++;
	      break;
	      
	    default:
	      Die("Traces2Alignment sees unrecognized statetype %c",
		  tr[seqidx]->statetype[statepos]);
	    }
	}
      for (; apos < alen; apos++) 
	aseqs[seqidx][apos] = '.';
      aseqs[seqidx][apos] = '\0';
    }

  /***********************************************
   * Build ainfo
   ***********************************************/
	
  if (ainfo != NULL)
    {
      ainfo->flags = AINFO_ALEN; 
      ainfo->alen  = alen;

      sprintf(ainfo->au, "HMM %s automatic alignment", RELEASE);
      ainfo->flags |= AINFO_AUTH;
      
			/* copy sqinfo structure array */
      if ((ainfo->sqinfo = (SQINFO *) malloc (sizeof(SQINFO) * nseq)) == NULL)
	Die("malloc failed");
      for (seqidx = 0; seqidx < nseq; seqidx++)
	SeqinfoCopy(&(ainfo->sqinfo[seqidx]), &(sqinfo[seqidx]));

      /* #=RF annotation:
       * x for match column, . for insert column
       */
      ainfo->flags |= AINFO_RF;
      ainfo->rf = (char *) MallocOrDie (sizeof(char) * (alen+1));
      apos = 0;
      if (!matchonly)
	for (count = 0; count < inserts[0]; count++)
	  ainfo->rf[apos++] = '.';
      else if (inserts[0] > 0)
	ainfo->rf[apos++] = '.';
      for (k = 1; k <= M; k++)
	{
	  ainfo->rf[apos++] = 'x';
	  if (!matchonly)
	    for (count = 0; count < inserts[k]; count++)
	      ainfo->rf[apos++] = '.';
	  else if (inserts[k] > 0)
	    ainfo->rf[apos++] = '.';
	}
      ainfo->rf[alen] = '\0';

      /* Currently, we produce no consensus structure. 
       * #=CS, generated from HMM internal annotation, would go here.
       */
    }

  /***********************************************
   * Garbage collection and return
   ***********************************************/
  free(inserts);
  *ret_aseqs = aseqs;
  return 1;
}
