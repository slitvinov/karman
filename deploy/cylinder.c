#line 1 "cylinder-cpp.c"
#line 1 "<built-in>"
#line 1 "<command-line>"
#line 1 "/usr/include/stdc-predef.h"
#line 1 "<command-line>"
#line 1 "cylinder-cpp.c"
#if _XOPEN_SOURCE < 700
  #undef _XOPEN_SOURCE
  #define _XOPEN_SOURCE 700
#endif
#if _GNU_SOURCE
#include <stdint.h>
#include <string.h>
#include <fenv.h>
#endif



#line 1 "/u/basilisk/src/common.h"
#line 1 "/u/basilisk/src/ast/std/stdlib.h"
#include <stdlib.h>
#line 2 "/u/basilisk/src/common.h"
#line 1 "/u/basilisk/src/ast/std/stdio.h"
#include <stdio.h>
#line 3 "/u/basilisk/src/common.h"
#line 1 "/u/basilisk/src/ast/std/stddef.h"
#include <stddef.h>
#line 4 "/u/basilisk/src/common.h"
#line 1 "/u/basilisk/src/ast/std/stdbool.h"
#include <stdbool.h>
#line 5 "/u/basilisk/src/common.h"
#line 1 "/u/basilisk/src/ast/std/stdarg.h"
#include <stdarg.h>
#line 6 "/u/basilisk/src/common.h"
#line 1 "/u/basilisk/src/ast/std/string.h"
#include <string.h>
#line 7 "/u/basilisk/src/common.h"
#line 1 "/u/basilisk/src/ast/std/float.h"
#include <float.h>
#line 8 "/u/basilisk/src/common.h"
#line 1 "/u/basilisk/src/ast/std/limits.h"
#include <limits.h>
#line 9 "/u/basilisk/src/common.h"
#line 1 "/u/basilisk/src/ast/std/math.h"
#include <math.h>
#line 10 "/u/basilisk/src/common.h"
#line 1 "/u/basilisk/src/ast/std/time.h"
#include <time.h>
#line 11 "/u/basilisk/src/common.h"
#line 1 "/u/basilisk/src/ast/std/sys/time.h"
#include <sys/time.h>
#line 12 "/u/basilisk/src/common.h"
#line 1 "/u/basilisk/src/ast/std/sys/resource.h"
#include <sys/resource.h>
#line 13 "/u/basilisk/src/common.h"

#if _OPENMP
# include <omp.h>
# define OMP(x) _Pragma(#x)
#elif _MPI

# define OMP(x)

# include <mpi.h>
static int mpi_rank, mpi_npe;
# define tid() mpi_rank
# define pid() mpi_rank
# define npe() mpi_npe

#else

# define OMP(x)

#endif
#line 46 "/u/basilisk/src/common.h"
#undef HUGE
#define HUGE ((double)1e30)

#define _NVARMAX 65536
#define is_constant(v) ((v).i >= _NVARMAX)
#define constant(v) (is_constant(v) ? _constant[(v).i - _NVARMAX] : HUGE)

#define max(a,b) ((a) > (b) ? (a) : (b))
#define min(a,b) ((a) < (b) ? (a) : (b))
#define sq(x) ((x)*(x))
#define cube(x) ((x)*(x)*(x))
#define sign(x) ((x) > 0 ? 1 : -1)
#define noise() (1. - 2.*rand()/(double)RAND_MAX)
#define clamp(x,a,b) ((x) < (a) ? (a) : (x) > (b) ? (b) : (x))

#define unmap(x,y)

#define trash(x)


#define systderr stderr
#define systdout stdout

#if _MPI
FILE * qstderr (void);
FILE * qstdout (void);
FILE * ferr = NULL, * fout = NULL;
#define not_mpi_compatible()\
do {\
  if (npe() > 1) {\
    fprintf (ferr, "%s() is not compatible with MPI (yet)\n", __func__);\
    exit (1);\
  }\
} while(0)\

#line 80

# define system(command) (pid() == 0 ? system(command) : 0)
#else
# define qstderr() stderr
# define qstdout() stdout
# define ferr stderr
# define fout stdout
# define not_mpi_compatible()
#endif



static inline void qassert (const char * file, int line, const char * cond) {
  fprintf (ferr, "%s:%d: Assertion `%s' failed.\n", file, line, cond);
  abort();
}
#line 104 "/u/basilisk/src/common.h"
#define sysmalloc malloc
#define syscalloc calloc
#define sysrealloc realloc
#define sysfree free
#define systrdup strdup

#if MTRACE

struct {
  FILE * fp;
  size_t total, max;
  size_t overhead, maxoverhead;
  size_t nr;
  size_t startrss, maxrss;
  char * fname;
} pmtrace;

typedef struct {
  char * func, * file;
  size_t max, total;
  int line, id;
} pmfunc;

typedef struct {
  size_t id, size;
} pmdata;

static pmfunc * pmfuncs = NULL;
static int pmfuncn = 0;

static int pmfunc_index (const char * func, const char * file, int line)
{
  pmfunc * p = pmfuncs;
  for (int i = 0; i < pmfuncn; i++, p++)
    if (p->line == line && !strcmp(func, p->func) && !strcmp(file, p->file))
      return p->id;
  pmfuncn++;
  pmfuncs = (pmfunc *) sysrealloc (pmfuncs, pmfuncn*sizeof(pmfunc));
  p = &pmfuncs[pmfuncn - 1];
  memset (p, 0, sizeof(pmfunc));
  p->func = systrdup(func);
  p->file = systrdup(file);
  p->line = line;
  p->id = pmfuncn;
  if (pmtrace.fp)
    fprintf (pmtrace.fp, "@ %d %s %s %d\n", pmfuncn, func, file, line);
  return pmfuncn;
}

static void pmfunc_trace (pmfunc * f, char c)
{
  if (pmtrace.fp)
    fprintf (pmtrace.fp, "%c %d %ld %ld %ld",
      c, f->id, pmtrace.nr, pmtrace.total, f->total);
#if _GNU_SOURCE
  if (pmtrace.nr % 1 == 0) {
    struct rusage usage;
    getrusage (RUSAGE_SELF, &usage);
    if (pmtrace.fp)
      fprintf (pmtrace.fp, " %ld", usage.ru_maxrss*1024);
    if (!pmtrace.nr)
      pmtrace.startrss = usage.ru_maxrss;
    if (usage.ru_maxrss > pmtrace.maxrss)
      pmtrace.maxrss = usage.ru_maxrss;
  }
#endif
  if (pmtrace.fp)
    fputc ('\n', pmtrace.fp);
  pmtrace.nr++;
}

static void * pmfunc_alloc (pmdata * d, size_t size,
       const char * func, const char * file, int line,
       char c)
{
  if (!(d != NULL)) qassert ("/u/basilisk/src/common.h", 0, "d != NULL");
  OMP (omp critical)
  {
    d->id = pmfunc_index(func, file, line);
    d->size = size;
    pmfunc * f = &pmfuncs[d->id - 1];
    f->total += size;
    if (f->total > f->max)
      f->max = f->total;
    pmtrace.total += size;
    pmtrace.overhead += sizeof(pmdata);
    if (pmtrace.total > pmtrace.max) {
      pmtrace.max = pmtrace.total;
      pmtrace.maxoverhead = pmtrace.overhead;
    }
    pmfunc_trace (f, c);
  }
  return ((char *)d) + sizeof(pmdata);
}

static void * pmfunc_free (void * ptr, char c)
{
  if (!ptr)
    return ptr;
  pmdata * d = (pmdata *) (((char *)ptr) - sizeof(pmdata));
  if (d->id < 1 || d->id > pmfuncn) {
    fputs ("*** MTRACE: ERROR!: corrupted free()", ferr);
    if (d->size == 0)
      fputs (", possible double free()", ferr);
    else
      fputs (", not traced?", ferr);
    fputs (", aborting...\n", ferr);
    abort();
    return ptr;
  }
  else
  OMP (omp critical)
  {
    pmfunc * f = &pmfuncs[d->id - 1];
    if (f->total < d->size) {
      fprintf (ferr, "*** MTRACE: ERROR!: %ld < %ld: corrupted free()?\n",
        f->total, d->size);
      abort();
    }
    else
      f->total -= d->size;
    if (pmtrace.total < d->size) {
      fprintf (ferr, "*** MTRACE: ERROR!: %ld < %ld: corrupted free()?\n",
        pmtrace.total, d->size);
      abort();
    }
    else {
      pmtrace.total -= d->size;
      pmtrace.overhead -= sizeof(pmdata);
    }
    d->id = 0;
    d->size = 0;
    pmfunc_trace (f, c);
  }
  return d;
}

static void * pmalloc (size_t size,
         const char * func, const char * file, int line)
{
  return pmfunc_alloc ((pmdata *) sysmalloc (sizeof(pmdata) + size),
         size, func, file, line, '+');
}

static void * pcalloc (size_t nmemb, size_t size,
         const char * func, const char * file, int line)
{
  void * p = pmalloc (nmemb*size, func, file, line);
  return memset (p, 0, nmemb*size);
}

static void * prealloc (void * ptr, size_t size,
   const char * func, const char * file, int line)
{
  return pmfunc_alloc ((pmdata *) sysrealloc (pmfunc_free(ptr, '<'),
           sizeof(pmdata) + size),
         size, func, file, line, '>');
}

static void pfree (void * ptr,
     const char * func, const char * file, int line)
{
  sysfree (pmfunc_free (ptr, '-'));
}

static char * pstrdup (const char * s,
         const char * func, const char * file, int line)
{
  char * d = (char *) pmalloc (strlen(s) + 1, func, file, line);
  return strcpy (d, s);
}

#if MTRACE < 3
static int pmaxsort (const void * a, const void * b) {
  const pmfunc * p1 = a, * p2 = b;
  return p1->max < p2->max;
}
#endif

static int ptotalsort (const void * a, const void * b) {
  const pmfunc * p1 = (const pmfunc *) a, * p2 = (const pmfunc *) b;
  return p1->total < p2->total;
}

static void pmfuncs_free()
{
  pmfunc * p = pmfuncs;
  for (int i = 0; i < pmfuncn; i++, p++) {
    sysfree (p->func);
    sysfree (p->file);
  }
  sysfree (pmfuncs);
}

void pmuntrace (void)
{
#if MTRACE < 3
  fprintf (ferr,
    "*** MTRACE: max resident  set size: %10ld bytes\n"
    "*** MTRACE: max traced memory size: %10ld bytes"
    " (tracing overhead %.1g%%)\n"
    "%10s    %20s   %s\n",
    pmtrace.maxrss*1024,
    pmtrace.max, pmtrace.maxoverhead*100./pmtrace.max,
    "max bytes", "function", "file");
  qsort (pmfuncs, pmfuncn, sizeof(pmfunc), pmaxsort);
  pmfunc * p = pmfuncs;
  for (int i = 0; i < pmfuncn && p->max > 0; i++, p++)
    fprintf (ferr, "%10ld    %20s   %s:%d\n",
      p->max, p->func, p->file, p->line);

  if (pmtrace.fp) {
    char * fname = pmtrace.fname, * s;
    while ((s = strchr(fname,'/')))
      fname = s + 1;

    fputs ("load(\"`echo $BASILISK`/mtrace.plot\")\n", pmtrace.fp);
    fprintf (pmtrace.fp,
      "plot '%s' u 3:($6-%g) w l t 'ru_maxrss - %.3g',"
      "total(\"%s\") w l t 'total'",
      fname,
      pmtrace.startrss*1024.,
      pmtrace.startrss*1024.,
      fname);
    pmfunc * p = pmfuncs;
    for (int i = 0; i < pmfuncn && p->max > 0.01*pmtrace.max; i++, p++)
      fprintf (pmtrace.fp,
        ",func(\"%s\",%d) w l t '%s'",
        fname, p->id, p->func);
    fputc ('\n', pmtrace.fp);
    fprintf (ferr,
      "*** MTRACE: To get a graph use: tail -n 2 %s | gnuplot -persist\n",
      fname);
    fclose (pmtrace.fp);
    pmtrace.fp = NULL;
    sysfree (pmtrace.fname);
  }
#endif

  if (pmtrace.total > 0) {
    qsort (pmfuncs, pmfuncn, sizeof(pmfunc), ptotalsort);
    pmfunc * p = pmfuncs;
    for (int i = 0; i < pmfuncn && p->total > 0; i++, p++)
      fprintf (ferr, "%s:%d: error: %ld bytes leaked here\n",
        p->file, p->line, p->total);
    pmfuncs_free();
    exit(1);
  }
  else {
#if MTRACE < 3
    fputs ("*** MTRACE: No memory leaks\n", ferr);
#endif
    pmfuncs_free();
  }
}

#else
# define pmalloc(s,func,file,line) malloc(s)
# define pcalloc(n,s,func,file,line) calloc(n,s)
# define prealloc(p,s,func,file,line) realloc(p,s)
# define pfree(p,func,file,line) free(p)
# define pstrdup(s,func,file,line) strdup(s)
#endif







typedef struct {
  void * p;
  long max, len;
} Array;

Array * array_new()
{
  Array * a = ((Array *) pmalloc ((1)*sizeof(Array),__func__,__FILE__,0));
  a->p = NULL;
  a->max = a->len = 0;
  return a;
}

void array_free (Array * a)
{
  pfree (a->p,__func__,__FILE__,0);
  pfree (a,__func__,__FILE__,0);
}

void array_append (Array * a, void * elem, size_t size)
{
  if (a->len + size >= a->max) {
    a->max += max (size, 4096);
    a->p = prealloc (a->p, a->max,__func__,__FILE__,0);
  }
  memcpy (((char *)a->p) + a->len, elem, size);
  a->len += size;
}

void * array_shrink (Array * a)
{
  void * p = prealloc (a->p, a->len,__func__,__FILE__,0);
  pfree (a,__func__,__FILE__,0);
  return p;
}



#if TRACE == 1
#include <extrae_user_events.h>

typedef struct {
  Array index, stack;
  extrae_type_t type;
} Trace;

Trace trace_func = {
  {NULL, 0, 0}, {NULL, 0, 0},
  60000010,
};

Trace trace_mpi_func = {
  {NULL, 0, 0}, {NULL, 0, 0},
  60000011,
};

static int lookup_func (Array * a, const char * func)
{
  for (int i = 0; i < a->len/sizeof(char *); i++) {
    char * s = ((char **)a->p)[i];
    if (!strcmp (func, s))
      return i + 1;
  }
  char * s = pstrdup (func,__func__,__FILE__,0);
  array_append (a, &s, sizeof(char *));
  return a->len;
}

static void trace_push (Trace * t, const char * func)
{
  int value = lookup_func (&t->index, func);
  Extrae_eventandcounters (t->type, value);
  array_append (&t->stack, &value, sizeof(int));
}

static void trace_pop (Trace * t, const char * func)
{
  if (!(t->stack.len > 0)) qassert ("/u/basilisk/src/common.h", 0, "t->stack.len > 0");
  t->stack.len -= sizeof(int);
  int value = t->stack.len > 0 ?
    ((int *)t->stack.p)[t->stack.len/sizeof(int) - 1] : 0;
  Extrae_eventandcounters (t->type, value);
}

static void trace_define (Trace * t, char * description)
{
  if (t->index.len > 0) {
    extrae_value_t values[t->index.len/sizeof(char *) + 1];
    char * names[t->index.len/sizeof(char *) + 1],
      ** func = (char **) t->index.p;
    names[0] = "OTHER";
    values[0] = 0;
    unsigned len = 1;
    for (int i = 0; i < t->index.len/sizeof(char *); i++, func++) {
      names[len] = *func;
      values[len++] = i + 1;
    }
    Extrae_define_event_type (&t->type, description, &len, values, names);
  }
}

static void trace_free (Trace * t)
{
  char ** func = (char **) t->index.p;
  for (int i = 0; i < t->index.len/sizeof(char *); i++, func++)
    pfree (*func,__func__,__FILE__,0);
  pfree (t->index.p,__func__,__FILE__,0);
  pfree (t->stack.p,__func__,__FILE__,0);
}

static void trace_off()
{
  trace_define (&trace_func, "Basilisk functions");
  trace_define (&trace_mpi_func, "Basilisk functions (MPI-related)");
  trace_free (&trace_func);
  trace_free (&trace_mpi_func);
}






# define tracing(func, file, line) trace_push (&trace_func, func)
# define end_tracing(func, file, line) trace_pop (&trace_func, func)

#elif TRACE

typedef struct {
  char * func, * file;
  int line, calls;
  double total, self;
#if _MPI
  double min, max;
#endif
} TraceIndex;

struct {
  Array stack, index;
  double t0;
} Trace = {
  {NULL, 0, 0}, {NULL, 0, 0},
  -1
};

static void trace_add (const char * func, const char * file, int line,
         double total, double self)
{
  TraceIndex * t = (TraceIndex *) Trace.index.p;
  int i, len = Trace.index.len/sizeof(TraceIndex);
  for (i = 0; i < len; i++, t++)
    if (t->line == line && !strcmp (func, t->func) && !strcmp (file, t->file))
      break;
  if (i == len) {
    TraceIndex t = {pstrdup(func,__func__,__FILE__,0), pstrdup(file,__func__,__FILE__,0), line, 1, total, self};
    array_append (&Trace.index, &t, sizeof(TraceIndex));
  }
  else
    t->calls++, t->total += total, t->self += self;
}

static void tracing (const char * func, const char * file, int line)
{
  struct timeval tv;
  gettimeofday (&tv, NULL);
  if (Trace.t0 < 0)
    Trace.t0 = tv.tv_sec + tv.tv_usec/1e6;
  double t[2] = {(tv.tv_sec - Trace.t0) + tv.tv_usec/1e6, 0.};
  array_append (&Trace.stack, t, 2*sizeof(double));




}

static void end_tracing (const char * func, const char * file, int line)
{
  struct timeval tv;
  gettimeofday (&tv, NULL);
  double te = (tv.tv_sec - Trace.t0) + tv.tv_usec/1e6;
  double * t = (double *) Trace.stack.p;
  if (!(Trace.stack.len >= 2*sizeof(double))) qassert ("/u/basilisk/src/common.h", 0, "Trace.stack.len >= 2*sizeof(double)");
  t += Trace.stack.len/sizeof(double) - 2;
  Trace.stack.len -= 2*sizeof(double);
  double dt = te - t[0];




  trace_add (func, file, line, dt, dt - t[1]);
  if (Trace.stack.len >= 2*sizeof(double)) {
    t -= 2;
    t[1] += dt;
  }
}

static int compar_self (const void * p1, const void * p2)
{
  const TraceIndex * t1 = p1, * t2 = p2;
  return t1->self < t2->self;
}

#if _MPI
static int compar_func (const void * p1, const void * p2)
{
  const TraceIndex * t1 = p1, * t2 = p2;
  if (t1->line != t2->line)
    return t1->line < t2->line;
  return strcmp (t1->file, t2->file);
}
#endif

void trace_print (FILE * fp, double threshold)
{
  int i, len = Trace.index.len/sizeof(TraceIndex);
  double total = 0.;
  TraceIndex * t;
  Array * index = array_new();
  for (i = 0, t = (TraceIndex *) Trace.index.p; i < len; i++, t++)
    array_append (index, t, sizeof(TraceIndex)), total += t->self;
#if _MPI
  qsort (index->p, len, sizeof(TraceIndex), compar_func);
  double tot[len], self[len], min[len], max[len];
  for (i = 0, t = (TraceIndex *) index->p; i < len; i++, t++)
    tot[i] = t->total, self[i] = t->self;
  MPI_Reduce (self, min, len, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
  MPI_Reduce (self, max, len, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce (pid() ? self : MPI_IN_PLACE,
       self, len, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce (pid() ? tot : MPI_IN_PLACE,
       tot, len, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  total = 0.;
  for (i = 0, t = (TraceIndex *) index->p; i < len; i++, t++)
    t->total = tot[i]/npe(), t->self = self[i]/npe(),
      t->max = max[i], t->min = min[i], total += t->self;
#endif
  qsort (index->p, len, sizeof(TraceIndex), compar_self);
  fprintf (fp, "   calls    total     self   %% total   function\n");
  for (i = 0, t = (TraceIndex *) index->p; i < len; i++, t++)
    if (t->self*100./total > threshold) {
      fprintf (fp, "%8d   %6.2f   %6.2f     %4.1f%%",
        t->calls, t->total, t->self, t->self*100./total);
#if _MPI
      fprintf (fp, " (%4.1f%% - %4.1f%%)", t->min*100./total, t->max*100./total);
#endif
      fprintf (fp, "   %s():%s:%d\n", t->func, t->file, t->line);
    }
  fflush (fp);
  array_free (index);
  for (i = 0, t = (TraceIndex *) Trace.index.p; i < len; i++, t++)
    t->calls = t->total = t->self = 0.;
}

static void trace_off()
{
  trace_print (fout, 0.);

  int i, len = Trace.index.len/sizeof(TraceIndex);
  TraceIndex * t;
  for (i = 0, t = (TraceIndex *) Trace.index.p; i < len; i++, t++)
    pfree (t->func,__func__,__FILE__,0), pfree (t->file,__func__,__FILE__,0);

  pfree (Trace.index.p,__func__,__FILE__,0);
  Trace.index.p = NULL;
  Trace.index.len = Trace.index.max = 0;

  pfree (Trace.stack.p,__func__,__FILE__,0);
  Trace.stack.p = NULL;
  Trace.stack.len = Trace.stack.max = 0;
}

#else
# define tracing(...)
# define end_tracing(...)
#endif



#if _OPENMP

#define tid() omp_get_thread_num()
#define pid() 0
#define npe() omp_get_num_threads()
#define mpi_all_reduce(v,type,op)
#define mpi_all_reduce_array(v,type,op,elem)

#elif _MPI

static bool in_prof = false;
static double prof_start, _prof;
#define prof_start(name)\
  if (!(!in_prof)) qassert ("/u/basilisk/src/common.h", 0, "!in_prof"); in_prof = true;\
  prof_start = MPI_Wtime();\

#line 667

#define prof_stop()\
  if (!(in_prof)) qassert ("/u/basilisk/src/common.h", 0, "in_prof"); in_prof = false;\
  _prof = MPI_Wtime();\
  mpi_time += _prof - prof_start;\

#line 672


#if FAKE_MPI
#define mpi_all_reduce(v,type,op)
#define mpi_all_reduce_array(v,type,op,elem)
#else
     
int mpi_all_reduce0 (void *sendbuf, void *recvbuf, int count,
       MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
{tracing("mpi_all_reduce0","/u/basilisk/src/common.h",0);
  { int _ret= MPI_Allreduce (sendbuf, recvbuf, count, datatype, op, comm);end_tracing("mpi_all_reduce0","/u/basilisk/src/common.h",0);return _ret;}
end_tracing("mpi_all_reduce0","/u/basilisk/src/common.h",0);}
#define mpi_all_reduce(v,type,op) {\
  prof_start ("mpi_all_reduce");\
  union { int a; float b; double c;} global;\
  mpi_all_reduce0 (&(v), &global, 1, type, op, MPI_COMM_WORLD);\
  memcpy (&(v), &global, sizeof (v));\
  prof_stop();\
}\

#line 691

#define mpi_all_reduce_array(v,type,op,elem) {\
  prof_start ("mpi_all_reduce");\
  type global[elem], tmp[elem];\
  for (int i = 0; i < elem; i++)\
    tmp[i] = (v)[i];\
  MPI_Datatype datatype;\
  if (!strcmp(#type, "double")) datatype = MPI_DOUBLE;\
  else if (!strcmp(#type, "int")) datatype = MPI_INT;\
  else if (!strcmp(#type, "long")) datatype = MPI_LONG;\
  else if (!strcmp(#type, "bool")) datatype = MPI_C_BOOL;\
  else {\
    fprintf (stderr, "unknown reduction type '%s'\n", #type);\
    fflush (stderr);\
    abort();\
  }\
  mpi_all_reduce0 (tmp, global, elem, datatype, op, MPI_COMM_WORLD);\
  for (int i = 0; i < elem; i++)\
    (v)[i] = global[i];\
  prof_stop();\
}\

#line 712


#endif

#define QFILE FILE

FILE * qstderr (void)
{
  static QFILE * fp = NULL;
  if (!fp) {
    if (mpi_rank > 0) {
      char name[80];
      sprintf (name, "log-%d", mpi_rank);
      fp = fopen (name, "w");
    }
    else
      fp = systderr;
  }
  return fp;
}

FILE * qstdout (void)
{
  static QFILE * fp = NULL;
  if (!fp) {
    if (mpi_rank > 0) {
      char name[80];
      sprintf (name, "out-%d", mpi_rank);
      fp = fopen (name, "w");
    }
    else
      fp = systdout;
  }
  return fp;
}

static void finalize (void)
{
  MPI_Finalize();
}

void mpi_init()
{
  int initialized;
  MPI_Initialized (&initialized);
  if (!initialized) {
    MPI_Init (NULL, NULL);
    MPI_Comm_set_errhandler (MPI_COMM_WORLD, MPI_ERRORS_ARE_FATAL);
    atexit (finalize);
  }
  MPI_Comm_rank (MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size (MPI_COMM_WORLD, &mpi_npe);
  srand (mpi_rank + 1);
  if (ferr == NULL) {
    if (mpi_rank > 0) {
      ferr = fopen ("/dev/null", "w");
      fout = fopen ("/dev/null", "w");
    }
    else {
      ferr = systderr;
      fout = systdout;
    }
    char * etrace = getenv ("MALLOC_TRACE"), name[80];
    if (etrace && mpi_rank > 0) {
      sprintf (name, "%s-%d", etrace, mpi_rank);
      setenv ("MALLOC_TRACE", name, 1);
    }
#if MTRACE == 1
    etrace = getenv ("MTRACE");
    if (!etrace)
      etrace = "mtrace";
    if (mpi_rank > 0) {
      sprintf (name, "%s-%d", etrace, mpi_rank);
      pmtrace.fp = fopen (name, "w");
      pmtrace.fname = systrdup(name);
    }
    else {
      pmtrace.fp = fopen (etrace, "w");
      pmtrace.fname = systrdup(etrace);
    }
#endif
  }
}

#else

#define tid() 0
#define pid() 0
#define npe() 1
#define mpi_all_reduce(v,type,op)
#define mpi_all_reduce_array(v,type,op,elem)

#endif

#define OMP_PARALLEL() OMP(omp parallel)

#define NOT_UNUSED(x) (void)(x)

#define VARIABLES ;
#define _index(a,m) (a.i)
#define val(a,k,l,m) data(k,l,m)[_index(a,m)]

double _val_higher_dimension = 0.;
#line 823 "/u/basilisk/src/common.h"
#if (_GNU_SOURCE || __APPLE__) && !_OPENMP && !_CADNA
double undefined;
# if __APPLE__
# include <stdint.h>
# include "fp_osx.h"
# endif
# define enable_fpe(flags) feenableexcept (flags)
# define disable_fpe(flags) fedisableexcept (flags)
static void set_fpe (void) {
  int64_t lnan = 0x7ff0000000000001;
  if (!(sizeof (int64_t) == sizeof (double))) qassert ("/u/basilisk/src/common.h", 0, "sizeof (int64_t) == sizeof (double)");
  memcpy (&undefined, &lnan, sizeof (double));
  enable_fpe (FE_DIVBYZERO|FE_INVALID);
}
#else
# define undefined ((double) DBL_MAX)
# define enable_fpe(flags)
# define disable_fpe(flags)
static void set_fpe (void) {}
#endif


typedef struct {
  long n;
  long tn;
  int depth;
  int maxdepth;
} Grid;
Grid * grid = NULL;

double X0 = 0., Y0 = 0., Z0 = 0.;

double L0 = 1.;


int N = 64;




typedef struct { int i; } scalar;

typedef struct {
  scalar x;

  scalar y;




} vector;

typedef struct {
  scalar * x;

  scalar * y;




} vectorl;

typedef struct {
  vector x;

  vector y;




} tensor;

struct { int x, y, z; } Period = {false, false, false};

typedef struct {
  double x, y, z;
} coord;

OMP(omp declare reduction (+ : coord :
      omp_out.x += omp_in.x,
      omp_out.y += omp_in.y,
      omp_out.z += omp_in.z))
#line 917 "/u/basilisk/src/common.h"
void normalize (coord * n)
{
  double norm = 0.;
  
    norm += sq(n->x);
    
#line 921
norm += sq(n->y);
  norm = sqrt(norm);
  
    n->x /= norm;
    
#line 924
n->y /= norm;
}

struct _origin { double x, y, z; };

void origin (struct _origin p) {
  X0 = p.x; Y0 = p.y; Z0 = p.z;
}

void size (double L) {
  L0 = L;
}

double zero (double s0, double s1, double s2) { return 0.; }






  enum { right, left, top, bottom };



int nboundary = 2*2;



#define dirichlet(expr) (2.*(expr) - val(_s,0,0,0))
#define dirichlet_homogeneous() (- val(_s,0,0,0))
#define dirichlet_face(expr) (expr)
#define dirichlet_face_homogeneous() (0.)
#define neumann(expr) (Delta*(expr) + val(_s,0,0,0))
#define neumann_homogeneous() (val(_s,0,0,0))

double * _constant = NULL;
size_t datasize = 0;
typedef struct _Point Point;

#line 1 "/u/basilisk/src/grid/boundaries.h"


typedef struct _Boundary Boundary;

struct _Boundary {
  void (* destroy) (Boundary * b);
  void (* level) (const Boundary * b, scalar * list, int l);

  void (* restriction) (const Boundary * b, scalar * list, int l);
};

static Boundary ** boundaries = NULL;

void add_boundary (Boundary * b) {
  int len = 0;
  if (boundaries) {
    Boundary ** i = boundaries;
    while (*i++) len++;
  }
  boundaries = (Boundary * *) prealloc (boundaries, (len + 2)*sizeof(Boundary *),__func__,__FILE__,0);
  boundaries[len] = b;
  boundaries[len+1] = NULL;
}

void free_boundaries() {
  if (!boundaries)
    return;
  Boundary ** i = boundaries, * b;
  while ((b = *i++))
    if (b->destroy)
      b->destroy (b);
    else
      pfree (b,__func__,__FILE__,0);
  pfree (boundaries,__func__,__FILE__,0);
  boundaries = NULL;
}
#line 47 "/u/basilisk/src/grid/boundaries.h"
typedef struct {
  Boundary parent;
  int d;
} BoxBoundary;
#line 964 "/u/basilisk/src/common.h"



typedef struct {
  double (** boundary) (Point, Point, scalar, void *);
  double (** boundary_homogeneous) (Point, Point, scalar, void *);
  double (* gradient) (double, double, double);
  void (* delete) (scalar);
  char * name;
  struct {
    int x;

    int y;




  } d;
  vector v;
  int face;
  bool nodump, freed;
  int block;
  scalar * depends;

  
#line 19 "/u/basilisk/src/grid/stencils.h"
bool input, output;
  int width;
  int dirty;
  
#line 18 "/u/basilisk/src/grid/multigrid-common.h"
void (* prolongation) (Point, scalar);
  void (* restriction) (Point, scalar);
  
#line 9 "/u/basilisk/src/grid/tree-common.h"
void (* refine) (Point, scalar);
  
#line 97
void (* coarsen) (Point, scalar);
  
#line 82 "/u/basilisk/src/fractions.h"
vector n;
  
#line 207 "/u/basilisk/src/embed-tree.h"
void (* embed_gradient) (Point, scalar, coord *);
  
#line 176 "/u/basilisk/src/embed.h"
bool third;

#line 987 "/u/basilisk/src/common.h"
} _Attributes;

static _Attributes * _attribute = NULL;






int list_len (scalar * list)
{
  if (!list) return 0;
  int ns = 0;
  {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ ns++;}}
  return ns;
}

scalar * list_append (scalar * list, scalar s)
{
  int len = list_len (list);
  list = (scalar *) prealloc (list, (len + 2)*sizeof(scalar),__func__,__FILE__,0);
  list[len] = s;
  list[len + 1].i = -1;
  return list;
}

scalar * list_prepend (scalar * list, scalar s)
{
  int len = list_len (list);
  list = (scalar *) prealloc (list, (len + 2)*sizeof(scalar),__func__,__FILE__,0);
  for (int i = len; i >= 1; i--)
    list[i] = list[i-1];
  list[0] = s;
  list[len + 1].i = -1;
  return list;
}

scalar * list_add (scalar * list, scalar s)
{
  {scalar*_i=(scalar*)( list);if(_i)for(scalar t=*_i;(&t)->i>=0;t=*++_i){
    if (t.i == s.i)
      return list;}}
  return list_append (list, s);
}

int list_lookup (scalar * l, scalar s)
{
  if (l != NULL)
    {scalar*_i=(scalar*)( l);if(_i)for(scalar s1=*_i;(&s1)->i>=0;s1=*++_i){
      if (s1.i == s.i)
 return true;}}
  return false;
}

scalar * list_copy (scalar * l)
{
  scalar * list = NULL;
  if (l != NULL)
    {scalar*_i=(scalar*)( l);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      list = list_append (list, s);}}
  return list;
}

scalar * list_concat (scalar * l1, scalar * l2)
{
  scalar * l3 = list_copy (l1);
  {scalar*_i=(scalar*)( l2);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    l3 = list_append (l3, s);}}
  return l3;
}

void list_print (scalar * l, FILE * fp)
{
  int i = 0;
  {scalar*_i=(scalar*)( l);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    fprintf (fp, "%s%s", i++ == 0 ? "{" : ",", _attribute[s.i].name);}}
  fputs (i > 0 ? "}\n" : "{}\n", fp);
}

int vectors_len (vector * list)
{
  if (!list) return 0;
  int nv = 0;
  {vector*_i=(vector*)( list);if(_i)for(vector v=*_i;(&v)->x.i>=0;v=*++_i){ nv++;}}
  return nv;
}

vector * vectors_append (vector * list, vector v)
{
  int len = vectors_len (list);
  list = (vector *) prealloc (list, (len + 2)*sizeof(vector),__func__,__FILE__,0);
  list[len] = v;
  list[len + 1] = (vector){{-1}};
  return list;
}

vector * vectors_add (vector * list, vector v)
{
  {vector*_i=(vector*)( list);if(_i)for(vector w=*_i;(&w)->x.i>=0;w=*++_i){ {
    bool id = true;
    
      if (w.x.i != v.x.i)
 id = false;
      
#line 1088
if (w.y.i != v.y.i)
 id = false;
    if (id)
      return list;
  }}}
  return vectors_append (list, v);
}

vector * vectors_copy (vector * l)
{
  vector * list = NULL;
  if (l != NULL)
    {vector*_i=(vector*)( l);if(_i)for(vector v=*_i;(&v)->x.i>=0;v=*++_i){
      list = vectors_append (list, v);}}
  return list;
}

vector * vectors_from_scalars (scalar * s)
{
  vector * list = NULL;
  while (s->i >= 0) {
    vector v;
     {
      if (!(s->i >= 0)) qassert ("/u/basilisk/src/common.h", 0, "s->i >= 0");
      v.x = *s++;
    } 
#line 1110
{
      if (!(s->i >= 0)) qassert ("/u/basilisk/src/common.h", 0, "s->i >= 0");
      v.y = *s++;
    }
    list = vectors_append (list, v);
  }
  return list;
}

int tensors_len (tensor * list)
{
  if (!list) return 0;
  int nt = 0;
  {tensor*_i=(tensor*)( list);if(_i)for(tensor t=*_i;(&t)->x.x.i>=0;t=*++_i){ nt++;}}
  return nt;
}

tensor * tensors_append (tensor * list, tensor t)
{
  int len = tensors_len (list);
  list = (tensor *) prealloc (list, (len + 2)*sizeof(tensor),__func__,__FILE__,0);
  list[len] = t;
  list[len + 1] = (tensor){{{-1}}};
  return list;
}

tensor * tensors_from_vectors (vector * v)
{
  tensor * list = NULL;
  while (v->x.i >= 0) {
    tensor t;
     {
      if (!(v->x.i >= 0)) qassert ("/u/basilisk/src/common.h", 0, "v->x.i >= 0");
      t.x = *v++;
    } 
#line 1141
{
      if (!(v->y.i >= 0)) qassert ("/u/basilisk/src/common.h", 0, "v->x.i >= 0");
      t.y = *v++;
    }
    list = tensors_append (list, t);
  }
  return list;
}

scalar * all = NULL;
scalar * baseblock = NULL;



scalar (* init_scalar) (scalar, const char *);
scalar (* init_vertex_scalar) (scalar, const char *);
vector (* init_vector) (vector, const char *);
tensor (* init_tensor) (tensor, const char *);
vector (* init_face_vector) (vector, const char *);





typedef struct _Event Event;
typedef int (* Expr) (int *, double *, Event *);

struct _Event {
  int last, nexpr;
  int (* action) (const int, const double, Event *);
  Expr expr[3];
  int * arrayi;
  double * arrayt;
  char * file;
  int line;
  char * name;
  double t;
  int i, a;
  void * data;
  Event * next;
};

static Event * Events = NULL;

int iter = 0, inext = 0;
double t = 0, tnext = 0;
void init_events (void);
void event_register (Event event);
static void _init_solver (void);

void init_solver()
{
  Events = pmalloc (sizeof (Event),__func__,__FILE__,0);
  Events[0].last = 1;
  _attribute = pcalloc (datasize/sizeof(double), sizeof (_Attributes),__func__,__FILE__,0);
  int n = datasize/sizeof(double);
  all = (scalar *) pmalloc (sizeof (scalar)*(n + 1),__func__,__FILE__,0);
  baseblock = (scalar *) pmalloc (sizeof (scalar)*(n + 1),__func__,__FILE__,0);
  for (int i = 0; i < n; i++)
    baseblock[i].i = all[i].i = i;
  baseblock[n].i = all[n].i = -1;
#if _CADNA
  cadna_init (-1);
#endif
#if _MPI
  mpi_init();
#elif MTRACE == 1
  char * etrace = getenv ("MTRACE");
  pmtrace.fp = fopen (etrace ? etrace : "mtrace", "w");
  pmtrace.fname = systrdup (etrace ? etrace : "mtrace");
#endif
}



#if _MPI
static double mpi_time = 0.;
#endif

typedef struct {
  clock_t c;
  struct timeval tv;
  double tm;
} timer;

timer timer_start (void)
{
  timer t;
  t.c = clock();
  gettimeofday (&t.tv, NULL);
#if _MPI
  t.tm = mpi_time;
#endif
  return t;
}

double timer_elapsed (timer t)
{
  struct timeval tvend;
  gettimeofday (&tvend, NULL);
  return ((tvend.tv_sec - t.tv.tv_sec) +
   (tvend.tv_usec - t.tv.tv_usec)/1e6);
}



const vector zerof = {{_NVARMAX+0},{_NVARMAX+1}};
const vector unityf = {{_NVARMAX+2},{_NVARMAX+3}};
const scalar unity = {_NVARMAX+4};
const scalar zeroc = {_NVARMAX+5};



        vector fm = {{_NVARMAX+2},{_NVARMAX+3}};
        scalar cm = {_NVARMAX+4};
#line 1269 "/u/basilisk/src/common.h"
static FILE ** qpopen_pipes = NULL;

FILE * qpopen (const char * command, const char * type)
{
  if (pid() > 0)
    return fopen ("/dev/null", type);
  FILE * fp = popen (command, type);
  if (fp) {
    FILE ** i = qpopen_pipes;
    int n = 0;
    while (i && *i) { n++; i++; }
    qpopen_pipes = (FILE * *) prealloc (qpopen_pipes, (n + 2)*sizeof(FILE *),__func__,__FILE__,0);
    qpopen_pipes[n] = fp;
    qpopen_pipes[n+1] = NULL;
  }
  return fp;
}

int qpclose (FILE * fp)
{
  if (pid() > 0)
    return fclose (fp);
  FILE ** i = qpopen_pipes;
  while (i && *i) {
    if (*i == fp)
      *i = (FILE *) 1;
    i++;
  }
  return pclose (fp);
}

static void qpclose_all()
{
  FILE ** i = qpopen_pipes;
  while (i && *i) {
    if (*i != (FILE *) 1)
      pclose (*i);
    i++;
  }
  pfree (qpopen_pipes,__func__,__FILE__,0);
  qpopen_pipes = NULL;
}






FILE * lfopen (const char * name, const char * mode)
{
  char fname[80];
  sprintf (fname, "%s-%d", name, pid());
  return fopen (fname, mode);
}



void * matrix_new (int n, int p, size_t size)
{
  void ** m = ((void * *) pmalloc ((n)*sizeof(void *),__func__,__FILE__,0));
  char * a = ((char *) pmalloc ((n*p*size)*sizeof(char),__func__,__FILE__,0));
  for (int i = 0; i < n; i++)
    m[i] = a + i*p*size;
  return m;
}

double matrix_inverse (double ** m, int n, double pivmin)
{
  int indxc[n], indxr[n], ipiv[n];
  int i, icol = 0, irow = 0, j, k, l, ll;
  double big, dum, pivinv, minpiv = HUGE;

  for (j = 0; j < n; j++)
    ipiv[j] = -1;

  for (i = 0; i < n; i++) {
    big = 0.0;
    for (j = 0; j < n; j++)
      if (ipiv[j] != 0)
 for (k = 0; k < n; k++) {
   if (ipiv[k] == -1) {
     if (fabs (m[j][k]) >= big) {
       big = fabs (m[j][k]);
       irow = j;
       icol = k;
     }
   }
 }
    ipiv[icol]++;
    if (irow != icol)
      for (l = 0; l < n; l++)
 do { double __tmp = m[irow][l]; m[irow][l] = m[icol][l]; m[icol][l] = __tmp; } while(0);
    indxr[i] = irow;
    indxc[i] = icol;
    if (fabs (m[icol][icol]) <= pivmin)
      return 0.;
    if (fabs (m[icol][icol]) < minpiv)
      minpiv = fabs (m[icol][icol]);
    pivinv = 1.0/m[icol][icol];
    m[icol][icol] = 1.0;
    for (l = 0; l < n; l++) m[icol][l] *= pivinv;
    for (ll = 0; ll < n; ll++)
      if (ll != icol) {
 dum = m[ll][icol];
 m[ll][icol] = 0.0;
 for (l = 0; l < n; l++)
   m[ll][l] -= m[icol][l]*dum;
      }
  }
  for (l = n - 1; l >= 0; l--) {
    if (indxr[l] != indxc[l])
      for (k = 0; k < n; k++)
 do { double __tmp = m[k][indxr[l]]; m[k][indxr[l]] = m[k][indxc[l]]; m[k][indxc[l]] = __tmp; } while(0);
  }
  return minpiv;
}

void matrix_free (void * m)
{
  pfree (((void **) m)[0],__func__,__FILE__,0);
  pfree (m,__func__,__FILE__,0);
}



typedef void (* free_solver_func) (void);

static Array * free_solver_funcs = NULL;

void free_solver_func_add (free_solver_func func)
{
  if (!free_solver_funcs)
    free_solver_funcs = array_new();
  array_append (free_solver_funcs, &func, sizeof(free_solver_func));
}



static char * display_defaults = NULL;

struct _display {
  const char * commands;
  bool overwrite;
};

static void free_display_defaults() {
  pfree (display_defaults,__func__,__FILE__,0);
}

void display (struct _display p)
{
  if (display_defaults == NULL)
    free_solver_func_add (free_display_defaults);
  if (p.overwrite) {
    pfree (display_defaults,__func__,__FILE__,0);
    display_defaults = pmalloc (strlen(p.commands) + 2,__func__,__FILE__,0);
    strcpy (display_defaults, "@");
    strcat (display_defaults, p.commands);
  }
  else {
    if (!display_defaults)
      display_defaults = pstrdup ("@",__func__,__FILE__,0);
    display_defaults =
      prealloc (display_defaults,
        strlen(display_defaults) + strlen(p.commands) + 1,__func__,__FILE__,0);
    strcat (display_defaults, p.commands);
  }
}







#line 1 "/u/basilisk/src/grid/stencils.h"
#line 17 "/u/basilisk/src/grid/stencils.h"










typedef struct {
  const char * fname;
  int line;
  int first;
  int face;
  bool vertex;
} ForeachData;


#define foreach_stencil() {\
  static ForeachData _loop = {\
    __FILE__, 0,\
    1, 0, 0\
  };\
  if (baseblock) for (scalar s = baseblock[0], * i = baseblock;\
  s.i >= 0; i++, s = *i) {\
    _attribute[s.i].input = _attribute[s.i].output = false;\
    _attribute[s.i].width = 0;\
  }\
  int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);\
  Point point = {0}; NOT_UNUSED (point);\

#line 48


#define end_foreach_stencil()\
  end_stencil (&_loop);\
  _loop.first = 0;\
}\

#line 54


#define foreach_vertex_stencil() foreach_stencil() _loop.vertex = true;
#define end_foreach_vertex_stencil() end_foreach_stencil()

#define foreach_face_stencil() foreach_stencil()
#define end_foreach_face_stencil() end_foreach_stencil()

#define foreach_visible_stencil(...) foreach_stencil()
#define end_foreach_visible_stencil(...) end_foreach_stencil()

#define _stencil_is_face_x() { _loop.face |= (1 << 0);
#define end__stencil_is_face_x() }
#define _stencil_is_face_y() { _loop.face |= (1 << 1);
#define end__stencil_is_face_y() }
#define _stencil_is_face_z() { _loop.face |= (1 << 2);
#define end__stencil_is_face_z() }

void stencil_val (Point p, scalar s, int i, int j, int k,
    const char * file, int line, bool overflow);
void stencil_val_a (Point p, scalar s, int i, int j, int k, bool input,
      const char * file, int line);

#define _stencil_val(a,_i,_j,_k)\
  stencil_val (point, a, _i, _j, _k, __FILE__, 0, false)\

#line 79

#define _stencil_val_o(a,_i,_j,_k)\
  stencil_val (point, a, _i, _j, _k, __FILE__, 0, true)\

#line 82

#define _stencil_val_a(a,_i,_j,_k)\
  stencil_val_a (point, a, _i, _j, _k, false, __FILE__, 0)\

#line 85

#define _stencil_val_r(a,_i,_j,_k)\
  stencil_val_a (point, a, _i, _j, _k, true, __FILE__, 0)\

#line 88


#define _stencil_fine(a,_i,_j,_k) _stencil_val(a,_i,_j,_k)
#define _stencil_fine_a(a,_i,_j,_k) _stencil_val_a(a,_i,_j,_k)
#define _stencil_fine_r(a,_i,_j,_k) _stencil_val_r(a,_i,_j,_k)

#define _stencil_coarse(a,_i,_j,_k) _stencil_val(a,_i,_j,_k)
#define _stencil_coarse_a(a,_i,_j,_k) _stencil_val_a(a,_i,_j,_k)
#define _stencil_coarse_r(a,_i,_j,_k) _stencil_val_r(a,_i,_j,_k)

#define r_assign(x)
#define _assign(x)

#define _stencil_neighbor(i,j,k)
#define _stencil_child(i,j,k)
#define _stencil_aparent(i,j,k)
#define _stencil_aparent_a(i,j,k)
#define _stencil_aparent_r(i,j,k)

#define _stencil_neighborp(i,j,k) neighborp(i,j,k)

int _stencil_nop;
#define _stencil_val_higher_dimension (_stencil_nop = 1)
#define _stencil__val_constant(a,_i,_j,_k) (_stencil_nop = 1)

typedef void _stencil_undefined;

#define o_stencil -2







static inline bool scalar_is_dirty (scalar s)
{
  if (_attribute[s.i].dirty)
    return true;
  scalar * depends = _attribute[s.i].depends;
  {scalar*_i=(scalar*)( depends);if(_i)for(scalar d=*_i;(&d)->i>=0;d=*++_i){
    if (_attribute[d.i].dirty)
      return true;}}
  return false;
}




static inline bool scalar_depends_from (scalar a, scalar b)
{
  scalar * depends = _attribute[a.i].depends;
  {scalar*_i=(scalar*)( depends);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    if (s.i == b.i)
      return true;}}
  return false;
}







void boundary_internal (scalar * list, const char * fname, int line);
void (* boundary_face) (vectorl);







void end_stencil (ForeachData * loop)
{
  scalar * listc = NULL, * dirty = NULL;
  vectorl listf = {NULL};
  bool flux = false;




  {scalar*_i=(scalar*)( baseblock);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {
    bool write = _attribute[s.i].output, read = _attribute[s.i].input;




    {





      if (read && scalar_is_dirty (s)) {





 if (_attribute[s.i].face) {
   if (_attribute[s.i].width > 0)
     listc = list_append (listc, s);
   else if (!write) {
     scalar sn = _attribute[s.i].v.x.i >= 0 ? _attribute[s.i].v.x : s;
     
       if (_attribute[s.i].v.x.i == s.i) {




  if (_attribute[sn.i].boundary[left] || _attribute[sn.i].boundary[right])
    listc = list_append (listc, s);
  else if (_attribute[s.i].dirty != 2) {
    listf.x = list_append (listf.x, s);
    flux = true;
  }
       }
       
#line 194
if (_attribute[s.i].v.y.i == s.i) {




  if (_attribute[sn.i].boundary[bottom] || _attribute[sn.i].boundary[top])
    listc = list_append (listc, s);
  else if (_attribute[s.i].dirty != 2) {
    listf.y = list_append (listf.y, s);
    flux = true;
  }
       }
   }
 }





 else if (_attribute[s.i].width > 0)
   listc = list_append (listc, s);
      }





      if (write) {
 if (2 > 1 && !loop->vertex && loop->first) {
   bool vertex = true;
   
     if (_attribute[s.i].d.x != -1)
       vertex = false;
     
#line 225
if (_attribute[s.i].d.y != -1)
       vertex = false;
   if (vertex)
     fprintf (ferr,
       "%s:%d: warning: vertex scalar '%s' should be assigned with"
       " a foreach_vertex() loop\n",
       loop->fname, loop->line, _attribute[s.i].name);
 }
 if (_attribute[s.i].face) {
   if (loop->face == 0 && loop->first)
     fprintf (ferr,
       "%s:%d: warning: face vector '%s' should be assigned with"
       " a foreach_face() loop\n",
       loop->fname, loop->line, _attribute[s.i].name);
 }
 else if (loop->face) {
   if (_attribute[s.i].v.x.i < 0) {
     int d = 1, i = 0;
      {
       if (loop->face == d) {
  _attribute[s.i].face = 2, _attribute[s.i].v.x.i = s.i;
  _attribute[s.i].boundary[left] = _attribute[s.i].boundary[right] = NULL;





       }
       d *= 2, i++;
     } 
#line 243
{
       if (loop->face == d) {
  _attribute[s.i].face = 2, _attribute[s.i].v.y.i = s.i;
  _attribute[s.i].boundary[bottom] = _attribute[s.i].boundary[top] = NULL;





       }
       d *= 2, i++;
     }
     if (!_attribute[s.i].face && loop->first)
       fprintf (ferr,
         "%s:%d: warning: scalar '%s' should be assigned with "
         "a foreach_face(x|y|z) loop\n",
         loop->fname, loop->line, _attribute[s.i].name);
   }
   else {
     char * name = NULL;
     if (_attribute[s.i].name) {
       name = pstrdup (_attribute[s.i].name,__func__,__FILE__,0);
       char * s = name + strlen(name) - 1;
       while (s != name && *s != '.') s--;
       if (s != name) *s = '\0';
     }
     struct { int x, y, z; } input, output;
     vector v = _attribute[s.i].v;

     
       input.x = _attribute[v.x.i].input, output.x = _attribute[v.x.i].output;
       
#line 273
input.y = _attribute[v.y.i].input, output.y = _attribute[v.y.i].output;

     init_face_vector (v, name);


     
       _attribute[v.x.i].input = input.x, _attribute[v.x.i].output = output.x;
       
#line 279
_attribute[v.y.i].input = input.y, _attribute[v.y.i].output = output.y;





     pfree (name,__func__,__FILE__,0);
   }
 }
 else if (loop->vertex) {
   bool vertex = true;
   
     if (_attribute[s.i].d.x != -1)
       vertex = false;
     
#line 291
if (_attribute[s.i].d.y != -1)
       vertex = false;
   if (!vertex) {
     char * name = NULL;
     if (_attribute[s.i].name) name = pstrdup (_attribute[s.i].name,__func__,__FILE__,0);
     init_vertex_scalar (s, name);
     
       _attribute[s.i].v.x.i = -1;
       
#line 298
_attribute[s.i].v.y.i = -1;




     pfree (name,__func__,__FILE__,0);
   }
 }





 dirty = list_append (dirty, s);
 {scalar*_i=(scalar*)( baseblock);if(_i)for(scalar d=*_i;(&d)->i>=0;d=*++_i){
   if (scalar_depends_from (d, s))
     dirty = list_append (dirty, d);}}
      }
    }
  }}}




  if (flux) {
#line 335 "/u/basilisk/src/grid/stencils.h"
    boundary_face (listf);
    
      pfree (listf.x,__func__,__FILE__,0);
      
#line 337
pfree (listf.y,__func__,__FILE__,0);
  }




  if (listc) {






    boundary_internal (listc, loop->fname, loop->line);
    pfree (listc,__func__,__FILE__,0);
  }





  if (dirty) {






    {scalar*_i=(scalar*)( dirty);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      _attribute[s.i].dirty = true;}}
    pfree (dirty,__func__,__FILE__,0);
  }
}
#line 1445 "/u/basilisk/src/common.h"
#line 14 "cylinder-cpp.c"
#line 1 "grid/quadtree.h"
#line 1 "/u/basilisk/src/grid/quadtree.h"


#line 1 "grid/tree.h"
#line 1 "/u/basilisk/src/grid/tree.h"
#line 1 "grid/mempool.h"
#line 1 "/u/basilisk/src/grid/mempool.h"





typedef struct _Pool Pool;

struct _Pool {
  Pool * next;
};

typedef struct {
  char * first, * lastb;
  size_t size;
  size_t poolsize;
  Pool * pool, * last;
} Mempool;

typedef struct {
  char * next;
} FreeBlock;

Mempool * mempool_new (size_t poolsize, size_t size)
{

  if (!(poolsize % 8 == 0)) qassert ("/u/basilisk/src/grid/mempool.h", 0, "poolsize % 8 == 0");
  if (!(size >= sizeof(FreeBlock))) qassert ("/u/basilisk/src/grid/mempool.h", 0, "size >= sizeof(FreeBlock)");


  poolsize = min(1 << 20, poolsize + sizeof(Pool));
  Mempool * m = ((Mempool *) pcalloc (1, sizeof(Mempool),__func__,__FILE__,0));
  m->poolsize = poolsize;
  m->size = size;
  return m;
}

void mempool_destroy (Mempool * m)
{
  Pool * p = m->pool;
  while (p) {
    Pool * next = p->next;
    pfree (p,__func__,__FILE__,0);
    p = next;
  }
  pfree (m,__func__,__FILE__,0);
}

void * mempool_alloc (Mempool * m)
{
  if (!m->first) {

    Pool * p = (Pool *) pmalloc (m->poolsize,__func__,__FILE__,0);
    p->next = NULL;
    if (m->last)
      m->last->next = p;
    else
      m->pool = p;
    m->last = p;
    m->first = m->lastb = ((char *)m->last) + sizeof(Pool);
    FreeBlock * b = (FreeBlock *) m->first;
    b->next = NULL;
  }
  void * ret = m->first;
  FreeBlock * b = (FreeBlock *) ret;
  char * next = b->next;
  if (!next) {
    m->lastb += m->size;
    next = m->lastb;
    if (next + m->size > ((char *) m->last) + m->poolsize)
      next = NULL;
    else {
      FreeBlock * b = (FreeBlock *) next;
      b->next = NULL;
    }
  }
  m->first = next;
#if TRASH
  double * v = (double *) ret;
  for (int i = 0; i < m->size/sizeof(double); i++)
    v[i] = undefined;
#endif
  return ret;
}

void * mempool_alloc0 (Mempool * m)
{
  void * ret = mempool_alloc (m);
  memset (ret, 0, m->size);
  return ret;
}

void mempool_free (Mempool * m, void * p)
{
#if TRASH
  double * v = (double *) p;
  for (int i = 0; i < m->size/sizeof(double); i++)
    v[i] = undefined;
#endif
  FreeBlock * b = (FreeBlock *) p;
  b->next = m->first;
  m->first = (char *) p;
}
#line 2 "/u/basilisk/src/grid/tree.h"




#line 1 "grid/memindex/range.h"
#line 1 "/u/basilisk/src/grid/memindex/range.h"
#line 15 "/u/basilisk/src/grid/memindex/range.h"
typedef struct {
  void ** p;
  int size;
} Memalloc;

typedef struct {
  int start, end;
} Memrange;
#line 34 "/u/basilisk/src/grid/memindex/range.h"
void memrange_alloc (Memrange * r, Memalloc * mem, int i)
{
  if (r->start == r->end) {
    r->start = i;
    r->end = i + 1;
    for (Memalloc * m = mem; m->p; m++) {
      *m->p = pcalloc (1, m->size,__func__,__FILE__,0);
      *m->p = (char *)(*m->p) - i*m->size;
    }
  }
  else if (i >= r->end) {
    for (Memalloc * m = mem; m->p; m++) {
      *m->p = prealloc ((char *)(*m->p) + r->start*m->size,
         m->size*(i + 1 - r->start),__func__,__FILE__,0);
      *m->p = (char *)(*m->p) - r->start*m->size;
      memset ((char *)(*m->p) + r->end*m->size, 0, (i - r->end + 1)*m->size);
    }
    r->end = i + 1;
  }
  else if (i < r->start) {
    for (Memalloc * m = mem; m->p; m++) {
      *m->p = prealloc ((char *)(*m->p) + r->start*m->size, m->size*(r->end - i),__func__,__FILE__,0);
      memmove ((char *)(*m->p) + (r->start - i)*m->size, *m->p,
        m->size*(r->end - r->start));
      memset ((char *)(*m->p), 0, (r->start - i)*m->size);
      *m->p = (char *)(*m->p) - i*m->size;
    }
    r->start = i;
  }
}
#line 73 "/u/basilisk/src/grid/memindex/range.h"
bool memrange_free (Memrange * r, Memalloc * mem, int i)
{
  if (i == r->start) {
    if (i == r->end - 1) {
      for (Memalloc * m = mem; m->p; m++) {
 pfree ((char *)(*m->p) + r->start*m->size,__func__,__FILE__,0);
 *m->p = NULL;
      }
      r->start = r->end = 0;
      return true;
    }
    else {
      for (i = i + 1; i < r->end &&
      !*(void **)((char *)(*mem->p) + i*mem->size); i++);
      for (Memalloc * m = mem; m->p; m++) {
 memmove ((char *)(*m->p) + r->start*m->size,
   (char *)(*m->p) + i*m->size, m->size*(r->end - i));
 *m->p = prealloc ((char *)(*m->p) + r->start*m->size,
    m->size*(r->end - i),__func__,__FILE__,0);
 *m->p = (char *)(*m->p) - i*m->size;
      }
      r->start = i;
    }
  }
  else if (i == r->end - 1) {
    for (i = i - 1; i >= r->start &&
    !*(void **)((char *)(*mem->p) + i*mem->size); i--);
    r->end = i + 1;
    for (Memalloc * m = mem; m->p; m++) {
      *m->p = prealloc ((char *)(*m->p) + r->start*m->size,
         m->size*(r->end - r->start),__func__,__FILE__,0);
      *m->p = (char *)(*m->p) - r->start*m->size;
    }
  }
  else {
    if (!(i > r->start && i < r->end)) qassert ("/u/basilisk/src/grid/memindex/range.h", 0, "i > r->start && i < r->end");
    for (Memalloc * m = mem; m->p; m++)
      memset ((char *)(*m->p) + i*m->size, 0, m->size);
  }
  return false;
}







struct _Memindex {
  Memrange r1;

  Memrange * r2;







  char *** b;



};
#line 171 "/u/basilisk/src/grid/memindex/range.h"
struct _Memindex * mem_new (int len)
{
  struct _Memindex * m = pcalloc (1, sizeof (struct _Memindex),__func__,__FILE__,0);
  return m;
}





void mem_destroy (struct _Memindex * m, int len)
{

  for (int i = m->r1.start; i < m->r1.end; i++)
    if (m->b[i]) {






      pfree (m->b[i] + m->r2[i].start,__func__,__FILE__,0);
    }
  if (m->b) {
    pfree (m->r2 + m->r1.start,__func__,__FILE__,0);



  }

  if (m->b)
    pfree (m->b + m->r1.start,__func__,__FILE__,0);
  pfree (m,__func__,__FILE__,0);
}
#line 218 "/u/basilisk/src/grid/memindex/range.h"
void mem_assign (struct _Memindex * m, int i, int j, int len, void * b)
{
  Memalloc mem[] = {{(void **)&m->b, sizeof(char **)},
      {(void **)&m->r2, sizeof(Memrange)},
      {NULL}};
  memrange_alloc (&m->r1, mem, i);
  Memalloc mem1[] = {{(void **)&m->b[i], sizeof(char *)},
       {NULL}};
  memrange_alloc (&m->r2[i], mem1, j);
  ((m)->b[i][j]) = b;
}
#line 259 "/u/basilisk/src/grid/memindex/range.h"
void mem_free (struct _Memindex * m, int i, int j, int len)
{
  Memalloc mem[] = {{(void **)&m->b[i], sizeof(char *)},
      {NULL}};
  if (memrange_free (&m->r2[i], mem, j)) {
    Memalloc mem[] = {{(void **)&m->b, sizeof(char **)},
        {(void **)&m->r2, sizeof(Memrange)},
        {NULL}};
    memrange_free (&m->r1, mem, i);
  }
}
#line 305 "/u/basilisk/src/grid/memindex/range.h"
#define foreach_mem(_m, _len, _i) {\
  Point point = {0};\
  for (point.i = max(Period.x*2, (_m)->r1.start);\
       point.i < min(_len - Period.x*2, (_m)->r1.end);\
       point.i += _i)\
    if ((_m)->b[point.i])\
      for (point.j = max(Period.y*2, (_m)->r2[point.i].start);\
    point.j < min(_len - Period.y*2, (_m)->r2[point.i].end);\
    point.j += _i)\
 if ((_m)->b[point.i][point.j]) {\

#line 315

#define end_foreach_mem() }}
#line 7 "/u/basilisk/src/grid/tree.h"
#line 24 "/u/basilisk/src/grid/tree.h"
typedef struct {
  unsigned short flags;

  unsigned short neighbors;
  int pid;
} Cell;

enum {
  active = 1 << 0,
  leaf = 1 << 1,
  border = 1 << 2,
  vertex = 1 << 3,
  user = 4,

  face_x = 1 << 0

  , face_y = 1 << 1




};

#define is_active(cell) ((cell).flags & active)
#define is_leaf(cell) ((cell).flags & leaf)
#define is_coarse() ((cell).neighbors > 0)
#define is_border(cell) ((cell).flags & border)
#define is_local(cell) ((cell).pid == pid())
#define is_vertex(cell) ((cell).flags & vertex)



typedef struct {
  int i;

  int j;




} IndexLevel;

typedef struct {
  IndexLevel * p;
  int n, nm;
} CacheLevel;

typedef struct {
  int i;

  int j;




  int level, flags;
} Index;

typedef struct {
  Index * p;
  int n, nm;
} Cache;



typedef struct {
  struct _Memindex * m;
  Mempool * pool;
  long nc;
  int len;
} Layer;

static size_t _size (size_t depth)
{
  return (1 << depth) + 2*2;
}

static size_t poolsize (size_t depth, size_t size)
{




  return sq(_size(depth))*size;



}

static Layer * new_layer (int depth)
{
  Layer * l = ((Layer *) pmalloc ((1)*sizeof(Layer),__func__,__FILE__,0));
  l->len = _size (depth);
  if (depth == 0)
    l->pool = NULL;
  else {
    size_t size = sizeof(Cell) + datasize;


    l->pool = mempool_new (poolsize (depth, size), (1 << 2)*size);
  }
  l->m = mem_new (l->len);
  l->nc = 0;
  return l;
}

static void destroy_layer (Layer * l)
{
  if (l->pool)
    mempool_destroy (l->pool);
  mem_destroy (l->m, l->len);
  pfree (l,__func__,__FILE__,0);
}



typedef struct {
  Grid g;
  Layer ** L;

  Cache leaves;
  Cache faces;
  Cache vertices;
  Cache refined;
  CacheLevel * active;
  CacheLevel * prolongation;
  CacheLevel * boundary;

  CacheLevel * restriction;

  bool dirty;
} Tree;



struct _Point {

  int i;

  int j;




  int level;
#ifdef foreach_block
  int l;
  #define _BLOCK_INDEX , point.l
#else
  #define _BLOCK_INDEX
#endif
};
static Point last_point;



static void cache_level_append (CacheLevel * c, Point p)
{
  if (c->n >= c->nm) {
    c->nm += 128;
    c->p = (IndexLevel *) prealloc (c->p, (c->nm)*sizeof(IndexLevel),__func__,__FILE__,0);
  }
  c->p[c->n].i = p.i;

  c->p[c->n].j = p.j;




  c->n++;
}

static void cache_level_shrink (CacheLevel * c)
{
  if (c->nm > (c->n/128 + 1)*128) {
    c->nm = (c->n/128 + 1)*128;
    if (!(c->nm > c->n)) qassert ("/u/basilisk/src/grid/tree.h", 0, "c->nm > c->n");
    c->p = (IndexLevel *) prealloc (c->p, sizeof (Index)*c->nm,__func__,__FILE__,0);
  }
}

static void cache_append (Cache * c, Point p, unsigned short flags)
{
  if (c->n >= c->nm) {
    c->nm += 128;
    c->p = (Index *) prealloc (c->p, (c->nm)*sizeof(Index),__func__,__FILE__,0);
  }
  c->p[c->n].i = p.i;

  c->p[c->n].j = p.j;




  c->p[c->n].level = p.level;
  c->p[c->n].flags = flags;
  c->n++;
}

void cache_shrink (Cache * c)
{
  cache_level_shrink ((CacheLevel *)c);
}
#line 243 "/u/basilisk/src/grid/tree.h"
#define allocated(k,l,n) (((point.i+k) >= (((Tree *)grid)->L[point.level]->m)->r1.start && (point.i+k) < (((Tree *)grid)->L[point.level]->m->r1.end) && (((Tree *)grid)->L[point.level]->m)->b[point.i+k] && (point.j+l) >= (((Tree *)grid)->L[point.level]->m)->r2[point.i+k].start && (point.j+l) < (((Tree *)grid)->L[point.level]->m)->r2[point.i+k].end && (((Tree *)grid)->L[point.level]->m)->b[point.i+k][point.j+l])\
                               )\

#line 245

#define NEIGHBOR(k,l,n) (((((Tree *)grid)->L[point.level]->m)->b[point.i+k][point.j+l])\
                            )\

#line 248

#define PARENT(k,l,n) (((((Tree *)grid)->L[point.level-1]->m)->b[(point.i+2)/2+k][(point.j+2)/2+l])\
                                                    )\

#line 251

#define allocated_child(k,l,n) (level < depth() &&\
         ((2*point.i-2 +k) >= (((Tree *)grid)->L[point.level+1]->m)->r1.start && (2*point.i-2 +k) < (((Tree *)grid)->L[point.level+1]->m->r1.end) && (((Tree *)grid)->L[point.level+1]->m)->b[2*point.i-2 +k] && (2*point.j-2 +l) >= (((Tree *)grid)->L[point.level+1]->m)->r2[2*point.i-2 +k].start && (2*point.j-2 +l) < (((Tree *)grid)->L[point.level+1]->m)->r2[2*point.i-2 +k].end && (((Tree *)grid)->L[point.level+1]->m)->b[2*point.i-2 +k][2*point.j-2 +l])\
\
                             )\

#line 256

#define CHILD(k,l,n) (((((Tree *)grid)->L[point.level+1]->m)->b[2*point.i-2 +k][2*point.j-2 +l])\
                                                )\

#line 259

#line 284 "/u/basilisk/src/grid/tree.h"
#define CELL(m) (*((Cell *)(m)))


#define depth() (grid->depth)
#define aparent(k,l,n) CELL(PARENT(k,l,n))
#define child(k,l,n) CELL(CHILD(k,l,n))


#define cell CELL(NEIGHBOR(0,0,0))
#define neighbor(k,l,n) CELL(NEIGHBOR(k,l,n))
#define neighborp(l,m,n) (Point) {\
    point.i + l,\
\
    point.j + m,\
\
\
\
\
    point.level\
    _BLOCK_INDEX\
}\

#line 305



#define data(k,l,n) ((double *) (NEIGHBOR(k,l,n) + sizeof(Cell)))
#define fine(a,k,p,n) ((double *) (CHILD(k,p,n) + sizeof(Cell)))[_index(a,n)]
#define coarse(a,k,p,n) ((double *) (PARENT(k,p,n) + sizeof(Cell)))[_index(a,n)]

#define POINT_VARIABLES\
  VARIABLES\
  int level = point.level; NOT_UNUSED(level);\
\
\
\
  struct { int x, y; } child = {\
    2*((point.i+2)%2)-1, 2*((point.j+2)%2)-1\
  };\
\
\
\
\
\
  NOT_UNUSED(child);\
  Point parent = point; NOT_UNUSED(parent);\
  parent.level--;\
  parent.i = (point.i + 2)/2;\
\
  parent.j = (point.j + 2)/2;\
\

#line 341


#line 1 "grid/foreach_cell.h"
#line 1 "/u/basilisk/src/grid/foreach_cell.h"
#line 66 "/u/basilisk/src/grid/foreach_cell.h"
#define foreach_cell_root(root)\
  {\
    int ig = 0, jg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg);\
    Point point = {0};\
\
\
\
    struct { int l, i, j, stage; } stack[20];\
\
\
\
\
    int _s = -1;\
    { _s++; stack[_s].l = 0; stack[_s].i = root.i; stack[_s].j = root.j; stack[_s].stage = 0; };\
    while (_s >= 0) {\
      int stage;\
      { point.level = stack[_s].l; point.i = stack[_s].i; point.j = stack[_s].j; stage = stack[_s].stage; _s--; };\
      if (!allocated (0,0,0))\
 continue;\
      switch (stage) {\
      case 0: {\
 POINT_VARIABLES;\
\

#line 89

#define end_foreach_cell_root()\
        if (point.level < grid->depth) {\
   { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].stage = 1; };\
          { _s++; stack[_s].l = point.level + 1; stack[_s].i = (2*point.i - 2); stack[_s].j = (2*point.j - 2); stack[_s].stage = 0; };\
        }\
        break;\
      }\
\
\
\
      case 1: { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].stage = 2; };\
       { _s++; stack[_s].l = point.level + 1; stack[_s].i = (2*point.i - 2); stack[_s].j = ((2*point.j - 2) + 1); stack[_s].stage = 0; }; break;\
      case 2: { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].stage = 3; };\
       { _s++; stack[_s].l = point.level + 1; stack[_s].i = ((2*point.i - 2) + 1); stack[_s].j = (2*point.j - 2); stack[_s].stage = 0; }; break;\
      case 3: { _s++; stack[_s].l = point.level + 1; stack[_s].i = ((2*point.i - 2) + 1); stack[_s].j = ((2*point.j - 2) + 1); stack[_s].stage = 0; }; break;\
\
      }\
    }\
  }\

#line 123


#define foreach_cell() {\
\
\
\
  Point root = {2,2,0};\
\
\
\
  foreach_cell_root (root)\

#line 134

#define end_foreach_cell() end_foreach_cell_root() }

#define foreach_cell_all() {\
  Point root = {0};\
  for (root.i = 2*Period.x; root.i <= 2*(2 - Period.x); root.i++)\
\
    for (root.j = 2*Period.y; root.j <= 2*(2 - Period.y); root.j++)\
\
\
\
\
 foreach_cell_root (root)\

#line 147

#define end_foreach_cell_all() end_foreach_cell_root() }

#define foreach_cell_post_root(condition, root)\
  {\
    int ig = 0, jg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg);\
    Point point = {0};\
\
\
\
    struct { int l, i, j, stage; } stack[20];\
\
\
\
\
    int _s = -1;\
    { _s++; stack[_s].l = 0; stack[_s].i = root.i; stack[_s].j = root.j; stack[_s].stage = 0; };\
    while (_s >= 0) {\
      int stage;\
      { point.level = stack[_s].l; point.i = stack[_s].i; point.j = stack[_s].j; stage = stack[_s].stage; _s--; };\
      if (!allocated (0,0,0))\
 continue;\
      switch (stage) {\
      case 0: {\
        POINT_VARIABLES;\
 if (point.level == grid->depth) {\
   { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].stage = 8; };\
 }\
 else {\
   { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].stage = 1; };\
   if (condition)\
     { _s++; stack[_s].l = point.level + 1; stack[_s].i = (2*point.i - 2); stack[_s].j = (2*point.j - 2); stack[_s].stage = 0; };\
 }\
 break;\
      }\
\
\
\
\
\
\
\
      case 1:\
 { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].stage = 2; };\
 if (condition)\
   { _s++; stack[_s].l = point.level + 1; stack[_s].i = (2*point.i - 2); stack[_s].j = ((2*point.j - 2) + 1); stack[_s].stage = 0; };\
 break;\
      case 2:\
 { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].stage = 3; };\
 if (condition)\
   { _s++; stack[_s].l = point.level + 1; stack[_s].i = ((2*point.i - 2) + 1); stack[_s].j = (2*point.j - 2); stack[_s].stage = 0; };\
 break;\
      case 3:\
 { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].stage = 4; };\
 if (condition)\
   { _s++; stack[_s].l = point.level + 1; stack[_s].i = ((2*point.i - 2) + 1); stack[_s].j = ((2*point.j - 2) + 1); stack[_s].stage = 0; };\
 break;\
\
      default: {\
        POINT_VARIABLES;\
\

#line 244

#define end_foreach_cell_post_root()\
      }\
      }\
    }\
  }\

#line 250


#define foreach_cell_post(condition)\
  {\
\
\
\
    Point root = {2,2,0};\
\
\
\
    foreach_cell_post_root(condition, root)\

#line 262

#define end_foreach_cell_post() end_foreach_cell_post_root() }

#define foreach_cell_post_all(condition) {\
  Point root = {0};\
  for (root.i = 0; root.i <= 2*2; root.i++)\
\
    for (root.j = 0; root.j <= 2*2; root.j++)\
\
\
\
\
 foreach_cell_post_root (condition, root)\

#line 275

#define end_foreach_cell_post_all() end_foreach_cell_post_root() }

#define foreach_leaf() foreach_cell()\
  if (is_leaf (cell)) {\
    if (is_active(cell) && is_local(cell)) {\

#line 281

#define end_foreach_leaf() } continue; } end_foreach_cell()
#line 344 "/u/basilisk/src/grid/tree.h"
#line 361 "/u/basilisk/src/grid/tree.h"
#define foreach_child() {\
  int _i = 2*point.i - 2, _j = 2*point.j - 2;\
  point.level++;\
  for (int _k = 0; _k < 2; _k++) {\
    point.i = _i + _k;\
    for (int _l = 0; _l < 2; _l++) {\
      point.j = _j + _l;\
      POINT_VARIABLES;\

#line 369

#define end_foreach_child()\
    }\
  }\
  point.i = (_i + 2)/2; point.j = (_j + 2)/2;\
  point.level--;\
}\

#line 376

#define foreach_child_break() _k = _l = 2
#line 407 "/u/basilisk/src/grid/tree.h"
#define is_refined_check() ((!is_leaf (cell) && cell.neighbors && cell.pid >= 0) &&\
    point.i > 0 && point.i < (1 << level) + 2*2 - 1\
\
    && point.j > 0 && point.j < (1 << level) + 2*2 - 1\
\
\
\
\
    )\

#line 416


#define foreach_cache(_cache) {\
  OMP_PARALLEL() {\
  int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);\
  Point point = {0};\
  point.i = 2;\
\
  point.j = 2;\
\
\
\
\
  int _k; unsigned short _flags; NOT_UNUSED(_flags);\
  OMP(omp for schedule(static))\
  for (_k = 0; _k < _cache.n; _k++) {\
    point.i = _cache.p[_k].i;\
\
    point.j = _cache.p[_k].j;\
\
\
\
\
    point.level = _cache.p[_k].level;\
    _flags = _cache.p[_k].flags;\
    POINT_VARIABLES;\

#line 442

#define end_foreach_cache() } } }

#define foreach_cache_level(_cache,_l) {\
  OMP_PARALLEL() {\
  int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);\
  Point point = {0};\
  point.i = 2;\
\
  point.j = 2;\
\
\
\
\
  point.level = _l;\
  int _k;\
  OMP(omp for schedule(static))\
  for (_k = 0; _k < _cache.n; _k++) {\
    point.i = _cache.p[_k].i;\
\
    point.j = _cache.p[_k].j;\
\
\
\
\
    POINT_VARIABLES;\

#line 468

#define end_foreach_cache_level() } } }

#define foreach_boundary_level(_l) {\
  if (_l <= depth()) {\
    { if (((Tree *)grid)->dirty) update_cache_f(); };\
    CacheLevel _boundary = ((Tree *)grid)->boundary[_l];\
    foreach_cache_level (_boundary,_l)\

#line 476

#define end_foreach_boundary_level() end_foreach_cache_level(); }}



#define foreach_boundary(_b) {\
  for (int _l = depth(); _l >= 0; _l--)\
    foreach_boundary_level(_l) {\
      if ((- cell.pid - 1) == _b)\
 for (int _d = 0; _d < 2; _d++) {\
   for (int _i = -1; _i <= 1; _i += 2) {\
     if (_d == 0) ig = _i; else if (_d == 1) jg = _i; else kg = _i;\
     if (allocated(-ig,-jg,-kg) &&\
  is_leaf (neighbor(-ig,-jg,-kg)) &&\
  !(neighbor(-ig,-jg,-kg).pid < 0) &&\
  is_local(neighbor(-ig,-jg,-kg))) {\
       point.i -= ig; x -= ig*Delta/2.;\
\
       point.j -= jg; y -= jg*Delta/2.;\
\
\
\
\

#line 499

#define end_foreach_boundary()\
       point.i += ig; x += ig*Delta/2.;\
\
       point.j += jg; y += jg*Delta/2.;\
\
\
\
\
            }\
   }\
   ig = jg = kg = 0;\
 }\
    } end_foreach_boundary_level(); }\

#line 513


#define foreach_halo(_name,_l) {\
  if (_l <= depth()) {\
    { if (((Tree *)grid)->dirty) update_cache_f(); };\
    CacheLevel _cache = ((Tree *)grid)->_name[_l];\
    foreach_cache_level (_cache, _l)\

#line 520

#define end_foreach_halo() end_foreach_cache_level(); }}

#line 1 "grid/neighbors.h"
#line 1 "/u/basilisk/src/grid/neighbors.h"
#line 17 "/u/basilisk/src/grid/neighbors.h"
#define foreach_neighbor(_s) {\
  int _nn = _s + 0 ? _s + 0 : 2;\
  int _i = point.i, _j = point.j;\
  for (int _k = - _nn; _k <= _nn; _k++) {\
    point.i = _i + _k;\
    for (int _l = - _nn; _l <= _nn; _l++) {\
      point.j = _j + _l;\
      POINT_VARIABLES;\

#line 25

#define end_foreach_neighbor()\
    }\
  }\
  point.i = _i; point.j = _j;\
}\

#line 31

#define foreach_neighbor_break() _k = _l = _nn + 1
#line 524 "/u/basilisk/src/grid/tree.h"

static inline bool has_local_children (Point point)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  {foreach_child()
    if (is_local(cell))
      return true;end_foreach_child()}
  return false;
}

static inline void cache_append_face (Point point, unsigned short flags)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  Tree * q = ((Tree *)grid);
  cache_append (&q->faces, point, flags);

  if (!is_vertex(cell)) {
    cache_append (&q->vertices, point, 0);
    cell.flags |= vertex;
  }
  
    if ((flags & face_y) && !is_vertex(neighbor(1,0,0))) {
      cache_append (&q->vertices, neighborp(1,0,0), 0);
      neighbor(1,0,0).flags |= vertex;
    }
    
#line 543
if ((flags & face_x) && !is_vertex(neighbor(0,1,0))) {
      cache_append (&q->vertices, neighborp(0,1,0), 0);
      neighbor(0,1,0).flags |= vertex;
    }
#line 557 "/u/basilisk/src/grid/tree.h"
}



static void update_cache_f (void)
{
  Tree * q = ((Tree *)grid);

  {foreach_cache (q->vertices)
    if (level <= depth() && allocated(0,0,0))
      cell.flags &= ~vertex;end_foreach_cache();}


  q->leaves.n = q->faces.n = q->vertices.n = 0;
  for (int l = 0; l <= depth(); l++)
    q->active[l].n = q->prolongation[l].n =
      q->boundary[l].n = q->restriction[l].n = 0;

  const unsigned short fboundary = 1 << user;
  {foreach_cell() {



    if (is_local(cell) && is_active(cell)) {


      cache_level_append (&q->active[level], point);
    }
#line 601 "/u/basilisk/src/grid/tree.h"
    if (!(cell.pid < 0)) {

      {foreach_neighbor (2)
 if (allocated(0,0,0) && (cell.pid < 0) && !(cell.flags & fboundary)) {
   cache_level_append (&q->boundary[level], point);
   cell.flags |= fboundary;
 }end_foreach_neighbor()}
    }

    else if (level > 0 && is_local(aparent(0,0,0)))
      cache_level_append (&q->restriction[level], point);

    if (is_leaf (cell)) {
      if (is_local(cell)) {
 cache_append (&q->leaves, point, 0);

 unsigned short flags = 0;
 
   if ((neighbor(-1,0,0).pid < 0) || (!is_leaf(neighbor(-1,0,0)) && !neighbor(-1,0,0).neighbors && neighbor(-1,0,0).pid >= 0) ||
       is_leaf(neighbor(-1,0,0)))
     flags |= face_x;
   
#line 619
if ((neighbor(0,-1,0).pid < 0) || (!is_leaf(neighbor(0,-1,0)) && !neighbor(0,-1,0).neighbors && neighbor(0,-1,0).pid >= 0) ||
       is_leaf(neighbor(0,-1,0)))
     flags |= face_y;
 if (flags)
   cache_append (&q->faces, point, flags);
 
   if ((neighbor(1,0,0).pid < 0) || (!is_leaf(neighbor(1,0,0)) && !neighbor(1,0,0).neighbors && neighbor(1,0,0).pid >= 0) ||
       (!is_local(neighbor(1,0,0)) && is_leaf(neighbor(1,0,0))))
     cache_append (&q->faces, neighborp(1,0,0), face_x);
   
#line 625
if ((neighbor(0,1,0).pid < 0) || (!is_leaf(neighbor(0,1,0)) && !neighbor(0,1,0).neighbors && neighbor(0,1,0).pid >= 0) ||
       (!is_local(neighbor(0,1,0)) && is_leaf(neighbor(0,1,0))))
     cache_append (&q->faces, neighborp(0,1,0), face_y);

 for (int i = 0; i <= 1; i++)

   for (int j = 0; j <= 1; j++)




       if (!is_vertex(neighbor(i,j,k))) {
  cache_append (&q->vertices, neighborp(i,j,k), 0);
  neighbor(i,j,k).flags |= vertex;
       }

        if (cell.neighbors > 0)
   cache_level_append (&q->prolongation[level], point);
      }
      else if (!(cell.pid < 0) || is_local(aparent(0,0,0))) {

 unsigned short flags = 0;
 
   if (allocated(-1,0,0) &&
       is_local(neighbor(-1,0,0)) && (!is_leaf(neighbor(-1,0,0)) && !neighbor(-1,0,0).neighbors && neighbor(-1,0,0).pid >= 0))
     flags |= face_x;
   
#line 648
if (allocated(0,-1,0) &&
       is_local(neighbor(0,-1,0)) && (!is_leaf(neighbor(0,-1,0)) && !neighbor(0,-1,0).neighbors && neighbor(0,-1,0).pid >= 0))
     flags |= face_y;
 if (flags)
   cache_append_face (point, flags);
 
   if (allocated(1,0,0) && is_local(neighbor(1,0,0)) &&
       (!is_leaf(neighbor(1,0,0)) && !neighbor(1,0,0).neighbors && neighbor(1,0,0).pid >= 0))
     cache_append_face (neighborp(1,0,0), face_x);
   
#line 654
if (allocated(0,1,0) && is_local(neighbor(0,1,0)) &&
       (!is_leaf(neighbor(0,1,0)) && !neighbor(0,1,0).neighbors && neighbor(0,1,0).pid >= 0))
     cache_append_face (neighborp(0,1,0), face_y);
      }

      continue;

    }
  }end_foreach_cell();}


  cache_shrink (&q->leaves);
  cache_shrink (&q->faces);
  cache_shrink (&q->vertices);
  for (int l = 0; l <= depth(); l++) {
    cache_level_shrink (&q->active[l]);
    cache_level_shrink (&q->prolongation[l]);
    cache_level_shrink (&q->boundary[l]);
    cache_level_shrink (&q->restriction[l]);
}

  q->dirty = false;


  for (int l = depth(); l >= 0; l--)
    {foreach_boundary_level (l)
      cell.flags &= ~fboundary;end_foreach_boundary_level();}



  grid->n = q->leaves.n;

#if !_MPI
  grid->tn = grid->n;
  grid->maxdepth = grid->depth;
#endif
}

#define foreach() { if (((Tree *)grid)->dirty) update_cache_f(); }; foreach_cache(((Tree *)grid)->leaves)
#define end_foreach() end_foreach_cache()

#define foreach_face_generic()\
  { if (((Tree *)grid)->dirty) update_cache_f(); };\
  foreach_cache(((Tree *)grid)->faces) 
#line 716

#define end_foreach_face_generic() end_foreach_cache()

#define is_face_x() { int ig = -1; VARIABLES; if (_flags & face_x) {
#define end_is_face_x() }}


#define is_face_y() { int jg = -1; VARIABLES; if (_flags & face_y) {
#define end_is_face_y() }}






#define foreach_vertex()\
  { if (((Tree *)grid)->dirty) update_cache_f(); };\
  foreach_cache(((Tree *)grid)->vertices) {\
    x -= Delta/2.;\
\
    y -= Delta/2.;\
\
\
\
\

#line 742

#define end_foreach_vertex() } end_foreach_cache()
#line 734 "/u/basilisk/src/grid/tree.h"
#define foreach_level(l) {\
  if (l <= depth()) {\
    { if (((Tree *)grid)->dirty) update_cache_f(); };\
    CacheLevel _active = ((Tree *)grid)->active[l];\
    foreach_cache_level (_active,l)\

#line 739

#define end_foreach_level() end_foreach_cache_level(); }}

#define foreach_coarse_level(l) foreach_level(l) if (!is_leaf(cell)) {
#define end_foreach_coarse_level() } end_foreach_level()

#define foreach_level_or_leaf(l) {\
  for (int _l1 = l; _l1 >= 0; _l1--)\
    foreach_level(_l1)\
      if (_l1 == l || is_leaf (cell)) {\

#line 749

#define end_foreach_level_or_leaf() } end_foreach_level(); }

#if TRASH
# undef trash
# define trash(list) reset(list, undefined)
#endif

void reset (void * alist, double val)
{
  scalar * list = (scalar *) alist;
  Tree * q = ((Tree *)grid);

  for (int l = 0; l <= depth(); l++) {
    Layer * L = q->L[l];
    {foreach_mem (L->m, L->len, 1) {
      point.level = l;
      {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {
 if (!is_constant(s))
   for (int b = 0; b < _attribute[s.i].block; b++)
     data(0,0,0)[s.i + b] = val;
      }}}
    }end_foreach_mem();}
  }
}

static CacheLevel * cache_level_resize (CacheLevel * name, int a)
{
  for (int i = 0; i <= depth() - a; i++)
    pfree (name[i].p,__func__,__FILE__,0);
  pfree (name,__func__,__FILE__,0);
  return ((CacheLevel *) pcalloc (depth() + 1, sizeof(CacheLevel),__func__,__FILE__,0));
}

static void update_depth (int inc)
{
  Tree * q = ((Tree *)grid);
  grid->depth += inc;
  q->L = &(q->L[-1]);
  q->L = (Layer * *) prealloc (q->L, (grid->depth + 2)*sizeof(Layer *),__func__,__FILE__,0);
  q->L = &(q->L[1]);
  if (inc > 0)
    q->L[grid->depth] = new_layer (grid->depth);
  q->active = cache_level_resize (q->active, inc);
  q->prolongation = cache_level_resize (q->prolongation, inc);
  q->boundary = cache_level_resize (q->boundary, inc);
  q->restriction = cache_level_resize (q->restriction, inc);
}
#line 823 "/u/basilisk/src/grid/tree.h"
typedef void (* PeriodicFunction) (struct _Memindex *, int, int, int, void *);

static void periodic_function (struct _Memindex * m, int i, int j, int len, void * b,
          PeriodicFunction f)
{
  f(m, i, j, len, b);
  if (Period.x) {
    int nl = len - 2*2;
    for (int l = - 1; l <= 1; l += 2)
      for (int n = i + l*nl; n >= 0 && n < len; n += l*nl)
 f(m, n, j, len, b);
    if (Period.y)
      for (int l = - 1; l <= 1; l += 2)
 for (int n = j + l*nl; n >= 0 && n < len; n += l*nl) {
   f(m, i, n, len, b);
   for (int o = - 1; o <= 1; o += 2)
     for (int p = i + o*nl; p >= 0 && p < len; p += o*nl)
       f(m, p, n, len, b);
 }
  }
  else if (Period.y) {
    int nl = len - 2*2;
    for (int l = - 1; l <= 1; l += 2)
      for (int n = j + l*nl; n >= 0 && n < len; n += l*nl)
 f(m, i, n, len, b);
  }
}

static void assign_periodic (struct _Memindex * m, int i, int j, int len, void * b)
{
  periodic_function (m, i, j, len, b, mem_assign);
}

static void free_periodic (struct _Memindex * m, int i, int j, int len)
{
  periodic_function (m, i, j, len, NULL, (PeriodicFunction) mem_free);
}
#line 938 "/u/basilisk/src/grid/tree.h"
static void alloc_children (Point point)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  if (point.level == grid->depth)
    update_depth (+1);
  else if (allocated_child(0,0,0))
    return;


  Layer * L = ((Tree *)grid)->L[point.level + 1];
  L->nc++;
  size_t len = sizeof(Cell) + datasize;
  char * b = (char *) mempool_alloc0 (L->pool);
  int i = 2*point.i - 2;
  for (int k = 0; k < 2; k++, i++) {




    int j = 2*point.j - 2;
    for (int l = 0; l < 2; l++, j++) {
      assign_periodic (L->m, i, j, L->len, b);
      b += len;
    }
#line 971 "/u/basilisk/src/grid/tree.h"
  }

  int pid = cell.pid;
  {foreach_child() {
    cell.pid = pid;
#if TRASH
    {scalar*_i=(scalar*)( all);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      val(s,0,0,0) = undefined;}}
#endif
  }end_foreach_child()}
}
#line 1000 "/u/basilisk/src/grid/tree.h"
static void free_children (Point point)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;

  Layer * L = ((Tree *)grid)->L[point.level + 1];
  int i = 2*point.i - 2, j = 2*point.j - 2;
  if (!(((L->m)->b[i][j]))) qassert ("/u/basilisk/src/grid/tree.h", 0, "mem_data (L->m,i,j)");
  mempool_free (L->pool, ((L->m)->b[i][j]));
  for (int k = 0; k < 2; k++)
    for (int l = 0; l < 2; l++)
      free_periodic (L->m, i + k, j + l, L->len);
  if (--L->nc == 0) {
    destroy_layer (L);
    if (!(point.level + 1 == grid->depth)) qassert ("/u/basilisk/src/grid/tree.h", 0, "point.level + 1 == grid->depth");
    update_depth (-1);
  }
}
#line 1041 "/u/basilisk/src/grid/tree.h"
void increment_neighbors (Point point)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  ((Tree *)grid)->dirty = true;
  if (cell.neighbors++ == 0)
    alloc_children (point);
  {foreach_neighbor (2/2)
    if (cell.neighbors++ == 0)
      alloc_children (point);end_foreach_neighbor()}
  cell.neighbors--;
}

void decrement_neighbors (Point point)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  ((Tree *)grid)->dirty = true;
  {foreach_neighbor (2/2)
    if (allocated(0,0,0)) {
      cell.neighbors--;
      if (cell.neighbors == 0)
 free_children (point);
    }end_foreach_neighbor()}
  if (cell.neighbors) {
    int pid = cell.pid;
    {foreach_child() {
      cell.flags = 0;
      cell.pid = pid;
    }end_foreach_child()}
  }
}

void realloc_scalar (int size)
{

  Tree * q = ((Tree *)grid);
  size_t oldlen = sizeof(Cell) + datasize;
  size_t newlen = oldlen + size;
  datasize += size;

  Layer * L = q->L[0];
  {foreach_mem (L->m, L->len, 1) {




    char * p = (char *) prealloc (((L->m)->b[point.i][point.j]),
     newlen*sizeof(char),__func__,__FILE__,0);
    assign_periodic (L->m, point.i, point.j, L->len, p);





  }end_foreach_mem();}

  for (int l = 1; l <= depth(); l++) {
    Layer * L = q->L[l];
    Mempool * oldpool = L->pool;
    L->pool = mempool_new (poolsize (l, newlen), (1 << 2)*newlen);
    {foreach_mem (L->m, L->len, 2) {
      char * new = (char *) mempool_alloc (L->pool);







      for (int k = 0; k < 2; k++)
 for (int o = 0; o < 2; o++) {
   memcpy (new, ((L->m)->b[point.i + k][point.j + o]), oldlen);
   assign_periodic (L->m, point.i + k, point.j + o, L->len, new);
   new += newlen;
 }
#line 1124 "/u/basilisk/src/grid/tree.h"
    }end_foreach_mem();}
    mempool_destroy (oldpool);
  }
}



#define VN v.x
#define VT v.y
#define VR v.z




#if _MPI
# define disable_fpe_for_mpi() disable_fpe (FE_DIVBYZERO|FE_INVALID)
# define enable_fpe_for_mpi() enable_fpe (FE_DIVBYZERO|FE_INVALID)
#else
# define disable_fpe_for_mpi()
# define enable_fpe_for_mpi()
#endif

static inline void no_restriction (Point point, scalar s);

static bool normal_neighbor (Point point, scalar * scalars, vector * vectors)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  for (int k = 1; k <= 2; k++)
    {
      for (int i = -k; i <= k; i += 2*k)
 if ((allocated(i,0,0) && !(neighbor(i,0,0).pid < 0))) {
   Point neighbor = neighborp(i,0,0);
   int id = (- cell.pid - 1);
   {scalar*_i=(scalar*)( scalars);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    
       val(s,0,0,0) = _attribute[s.i].boundary[id](neighbor, point, s, NULL);}}
   {vector*_i=(vector*)( vectors);if(_i)for(vector v=*_i;(&v)->x.i>=0;v=*++_i){
     {
       scalar vn = VN;
       val(v.x,0,0,0) = _attribute[vn.i].boundary[id](neighbor, point, v.x, NULL);

       scalar vt = VT;
       val(v.y,0,0,0) = _attribute[vt.i].boundary[id](neighbor, point, v.y, NULL);





     }}}
   return true;
 }
      
#line 1152
for (int i = -k; i <= k; i += 2*k)
 if ((allocated(0,i,0) && !(neighbor(0,i,0).pid < 0))) {
   Point neighbor = neighborp(0,i,0);
   int id = (- cell.pid - 1);
   {scalar*_i=(scalar*)( scalars);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    
       val(s,0,0,0) = _attribute[s.i].boundary[id](neighbor, point, s, NULL);}}
   {vector*_i=(vector*)( vectors);if(_i)for(vector v=*_i;(&v)->x.i>=0;v=*++_i){
     {
       scalar vn = VN;
       val(v.y,0,0,0) = _attribute[vn.i].boundary[id](neighbor, point, v.y, NULL);

       scalar vt = VT;
       val(v.x,0,0,0) = _attribute[vt.i].boundary[id](neighbor, point, v.x, NULL);





     }}}
   return true;
 }}
  return false;
}

static bool diagonal_neighbor_2D (Point point,
      scalar * scalars, vector * vectors)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;

  for (int k = 1; k <= 2; k++)



      for (int i = -k; i <= k; i += 2*k)
 for (int j = -k; j <= k; j += 2*k)
   if (allocated(i,j,0) && (allocated(i,j,0) && !(neighbor(i,j,0).pid < 0)) &&
       allocated(i,0,0) && (neighbor(i,0,0).pid < 0) &&
       allocated(0,j,0) && (neighbor(0,j,0).pid < 0)) {
     Point n = neighborp(i,j,0),
       n1 = neighborp(i,0,0), n2 = neighborp(0,j,0);
     int id1 = (- neighbor(i,0,0).pid - 1), id2 = (- neighbor(0,j,0).pid - 1);
     {scalar*_i=(scalar*)( scalars);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      
  val(s,0,0,0) = (_attribute[s.i].boundary[id1](n,n1,s,NULL) +
         _attribute[s.i].boundary[id2](n,n2,s,NULL) -
         val(s,i,j,0));}}
     {vector*_i=(vector*)( vectors);if(_i)for(vector v=*_i;(&v)->x.i>=0;v=*++_i){
       {
  scalar vt = VT, vn = VN;
  val(v.x,0,0,0) = (_attribute[vt.i].boundary[id1](n,n1,v.x,NULL) +
    _attribute[vn.i].boundary[id2](n,n2,v.x,NULL) -
    val(v.x,i,j,0));
  val(v.y,0,0,0) = (_attribute[vn.i].boundary[id1](n,n1,v.y,NULL) +
    _attribute[vt.i].boundary[id2](n,n2,v.y,NULL) -
    val(v.y,i,j,0));






       }}}
     return true;
   }

  return false;
}

static bool diagonal_neighbor_3D (Point point,
      scalar * scalars, vector * vectors)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
#line 1266 "/u/basilisk/src/grid/tree.h"
  return false;
}



static Point tangential_neighbor_x (Point point, bool * zn)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  for (int k = 1; k <= 2; k++)
    for (int j = -k; j <= k; j += 2*k) {
      if ((allocated(0,j,0) && !(neighbor(0,j,0).pid < 0)) || (allocated(-1,j,0) && !(neighbor(-1,j,0).pid < 0))) {
 *zn = false;
 return neighborp(0,j,0);
      }







    }
  return (Point){.level = -1};
}

#line 1271
static Point tangential_neighbor_y (Point point, bool * zn)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  for (int k = 1; k <= 2; k++)
    for (int j = -k; j <= k; j += 2*k) {
      if ((allocated(j,0,0) && !(neighbor(j,0,0).pid < 0)) || (allocated(j,-1,0) && !(neighbor(j,-1,0).pid < 0))) {
 *zn = false;
 return neighborp(j,0,0);
      }







    }
  return (Point){.level = -1};
}


static inline bool is_boundary_point (Point point) {int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  return (cell.pid < 0);
}

static void box_boundary_level (const Boundary * b, scalar * list, int l)
{
  disable_fpe_for_mpi();
  scalar * scalars = NULL;
  vector * vectors = NULL, * faces = NULL;
  {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    if (!is_constant(s) && _attribute[s.i].refine != no_restriction) {
      if (_attribute[s.i].v.x.i == s.i) {
 if (_attribute[s.i].face)
   faces = vectors_add (faces, _attribute[s.i].v);
 else
   vectors = vectors_add (vectors, _attribute[s.i].v);
      }
      else if (_attribute[s.i].v.x.i < 0 && _attribute[s.i].boundary[0])
 scalars = list_add (scalars, s);
    }}}

  {foreach_boundary_level (l) {
    if (!normal_neighbor (point, scalars, vectors) &&
 !diagonal_neighbor_2D (point, scalars, vectors) &&
 !diagonal_neighbor_3D (point, scalars, vectors)) {

      {scalar*_i=(scalar*)( scalars);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){

   val(s,0,0,0) = undefined;}}
      {vector*_i=(vector*)( vectors);if(_i)for(vector v=*_i;(&v)->x.i>=0;v=*++_i){

   {
     val(v.x,0,0,0) = undefined;
     
#line 1323
val(v.y,0,0,0) = undefined;}}}
    }
    if (faces) {
      int id = (- cell.pid - 1);
      
 for (int i = -1; i <= 1; i += 2) {

   if ((allocated(i,0,0) && !(neighbor(i,0,0).pid < 0))) {
     Point neighbor = neighborp(i,0,0);
     {vector*_i=(vector*)( faces);if(_i)for(vector v=*_i;(&v)->x.i>=0;v=*++_i){ {
       scalar vn = VN;
       if (_attribute[vn.i].boundary[id])
 
    val(v.x,(i + 1)/2,0,0) = _attribute[vn.i].boundary[id](neighbor, point, v.x, NULL);
     }}}
   }

   else if (i == -1) {

     bool zn;
     Point neighbor = tangential_neighbor_x (point, &zn);
     if (neighbor.level >= 0) {
       int id = is_boundary_point (neighbor) ?
  (- neighbor(-1,0,0).pid - 1) : (- cell.pid - 1);
       {vector*_i=(vector*)( faces);if(_i)for(vector v=*_i;(&v)->x.i>=0;v=*++_i){ {

  scalar vt = VT;



 
    val(v.x,0,0,0) = _attribute[vt.i].boundary[id](neighbor, point, v.x, NULL);
       }}}
     }
     else

       {vector*_i=(vector*)( faces);if(_i)for(vector v=*_i;(&v)->x.i>=0;v=*++_i){
 
    val(v.x,0,0,0) = 0.;}}
   }

 }
 
#line 1328
for (int i = -1; i <= 1; i += 2) {

   if ((allocated(0,i,0) && !(neighbor(0,i,0).pid < 0))) {
     Point neighbor = neighborp(0,i,0);
     {vector*_i=(vector*)( faces);if(_i)for(vector v=*_i;(&v)->x.i>=0;v=*++_i){ {
       scalar vn = VN;
       if (_attribute[vn.i].boundary[id])
 
    val(v.y,0,(i + 1)/2,0) = _attribute[vn.i].boundary[id](neighbor, point, v.y, NULL);
     }}}
   }

   else if (i == -1) {

     bool zn;
     Point neighbor = tangential_neighbor_y (point, &zn);
     if (neighbor.level >= 0) {
       int id = is_boundary_point (neighbor) ?
  (- neighbor(0,-1,0).pid - 1) : (- cell.pid - 1);
       {vector*_i=(vector*)( faces);if(_i)for(vector v=*_i;(&v)->x.i>=0;v=*++_i){ {

  scalar vt = VT;



 
    val(v.y,0,0,0) = _attribute[vt.i].boundary[id](neighbor, point, v.y, NULL);
       }}}
     }
     else

       {vector*_i=(vector*)( faces);if(_i)for(vector v=*_i;(&v)->x.i>=0;v=*++_i){
 
    val(v.y,0,0,0) = 0.;}}
   }

 }
    }
  }end_foreach_boundary_level();}

  pfree (scalars,__func__,__FILE__,0);
  pfree (vectors,__func__,__FILE__,0);
  pfree (faces,__func__,__FILE__,0);
  enable_fpe_for_mpi();
}



#undef VN
#undef VT
#define VN _attribute[s.i].v.x
#define VT _attribute[s.i].v.y

static double masked_average (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  double sum = 0., n = 0.;
  {foreach_child()
    if (!(cell.pid < 0) && val(s,0,0,0) != HUGE)
      sum += val(s,0,0,0), n++;end_foreach_child()}
  return n ? sum/n : HUGE;
}


static double masked_average_x (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  double sum = 0., n = 0.;
  {foreach_child()
    if (child.x < 0 && (!(cell.pid < 0) || !(neighbor(1,0,0).pid < 0)) &&
 val(s,1,0,0) != HUGE)
      sum += val(s,1,0,0), n++;end_foreach_child()}
  return n ? sum/n : HUGE;
}

#line 1391
static double masked_average_y (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  double sum = 0., n = 0.;
  {foreach_child()
    if (child.y < 0 && (!(cell.pid < 0) || !(neighbor(0,1,0).pid < 0)) &&
 val(s,0,1,0) != HUGE)
      sum += val(s,0,1,0), n++;end_foreach_child()}
  return n ? sum/n : HUGE;
}

static void masked_boundary_restriction (const Boundary * b,
      scalar * list, int l)
{
  scalar * scalars = NULL;
  vector * faces = NULL;
  {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    if (!is_constant(s) && _attribute[s.i].refine != no_restriction) {
      if (_attribute[s.i].v.x.i == s.i && _attribute[s.i].face)
 faces = vectors_add (faces, _attribute[s.i].v);
      else
 scalars = list_add (scalars, s);
    }}}

  {foreach_halo (restriction, l) {
    {scalar*_i=(scalar*)( scalars);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      val(s,0,0,0) = masked_average (parent, s);}}
    {vector*_i=(vector*)( faces);if(_i)for(vector v=*_i;(&v)->x.i>=0;v=*++_i){
      { {
 double average = masked_average_x (parent, v.x);
 if ((neighbor(-1,0,0).pid < 0))
   val(v.x,0,0,0) = average;
 if ((neighbor(1,0,0).pid < 0))
   val(v.x,1,0,0) = average;
      } 
#line 1418
{
 double average = masked_average_y (parent, v.y);
 if ((neighbor(0,-1,0).pid < 0))
   val(v.y,0,0,0) = average;
 if ((neighbor(0,1,0).pid < 0))
   val(v.y,0,1,0) = average;
      }}}}
  }end_foreach_halo();}

  pfree (scalars,__func__,__FILE__,0);
  pfree (faces,__func__,__FILE__,0);
}
#line 1454 "/u/basilisk/src/grid/tree.h"
static void free_cache (CacheLevel * c)
{
  for (int l = 0; l <= depth(); l++)
    pfree (c[l].p,__func__,__FILE__,0);
  pfree (c,__func__,__FILE__,0);
}

void free_grid (void)
{
  if (!grid)
    return;
  free_boundaries();
  Tree * q = ((Tree *)grid);
  pfree (q->leaves.p,__func__,__FILE__,0);
  pfree (q->faces.p,__func__,__FILE__,0);
  pfree (q->vertices.p,__func__,__FILE__,0);
  pfree (q->refined.p,__func__,__FILE__,0);


  Layer * L = q->L[0];
  {foreach_mem (L->m, L->len, 1) {



    pfree (((L->m)->b[point.i][point.j]),__func__,__FILE__,0);



  }end_foreach_mem();}
  for (int l = 0; l <= depth(); l++)
    destroy_layer (q->L[l]);
  q->L = &(q->L[-1]);
  pfree (q->L,__func__,__FILE__,0);
  free_cache (q->active);
  free_cache (q->prolongation);
  free_cache (q->boundary);
  free_cache (q->restriction);
  pfree (q,__func__,__FILE__,0);
  grid = NULL;
}

static void refine_level (int depth);

     
void init_grid (int n)
{tracing("init_grid","/u/basilisk/src/grid/tree.h",0);

  if (!(sizeof(Cell) % 8 == 0)) qassert ("/u/basilisk/src/grid/tree.h", 0, "sizeof(Cell) % 8 == 0");

  free_grid();
  int depth = 0;
  while (n > 1) {
    if (n % 2) {
      fprintf (ferr, "tree: N must be a power-of-two\n");
      exit (1);
    }
    n /= 2;
    depth++;
  }
  Tree * q = ((Tree *) pcalloc (1, sizeof(Tree),__func__,__FILE__,0));
  grid = (Grid *) q;
  grid->depth = 0;


  q->L = ((Layer * *) pmalloc ((2)*sizeof(Layer *),__func__,__FILE__,0));

  q->L[0] = NULL; q->L = &(q->L[1]);

  Layer * L = new_layer (0);
  q->L[0] = L;
#line 1537 "/u/basilisk/src/grid/tree.h"
  for (int i = Period.x*2; i < L->len - Period.x*2; i++)
    for (int j = Period.y*2; j < L->len - Period.y*2; j++)
      assign_periodic (L->m, i, j, L->len,
         (char *) pcalloc (1, sizeof(Cell) + datasize,__func__,__FILE__,0));
  CELL(((L->m)->b[2][2])).flags |= leaf;
  if (pid() == 0)
    CELL(((L->m)->b[2][2])).flags |= active;
  for (int k = - 2*(1 - Period.x); k <= 2*(1 - Period.x); k++)
    for (int l = -2*(1 - Period.y); l <= 2*(1 - Period.y); l++)
      CELL(((L->m)->b[2 +k][2 +l])).pid =
 (k < 0 ? -1 - left :
  k > 0 ? -1 - right :
  l > 0 ? -1 - top :
  l < 0 ? -1 - bottom :
  0);
  CELL(((L->m)->b[2][2])).pid = 0;
#line 1575 "/u/basilisk/src/grid/tree.h"
  q->active = ((CacheLevel *) pcalloc (1, sizeof(CacheLevel),__func__,__FILE__,0));
  q->prolongation = ((CacheLevel *) pcalloc (1, sizeof(CacheLevel),__func__,__FILE__,0));
  q->boundary = ((CacheLevel *) pcalloc (1, sizeof(CacheLevel),__func__,__FILE__,0));
  q->restriction = ((CacheLevel *) pcalloc (1, sizeof(CacheLevel),__func__,__FILE__,0));
  q->dirty = true;
  N = 1 << depth;
#if _MPI
  void mpi_boundary_new();
  mpi_boundary_new();
#endif

  Boundary * b = ((Boundary *) pcalloc (1, sizeof(Boundary),__func__,__FILE__,0));
  b->level = box_boundary_level;
  b->restriction = masked_boundary_restriction;
  add_boundary (b);
  refine_level (depth);
  reset (all, 0.);
  { if (((Tree *)grid)->dirty) update_cache_f(); };
end_tracing("init_grid","/u/basilisk/src/grid/tree.h",0);}


void check_two_one (void)
{
  {foreach_leaf()
    if (level > 0)
      for (int k = -1; k <= 1; k++)
 for (int l = -1; l <= 1; l++) {

   int i = (point.i + 2)/2 + k;
   int j = (point.j + 2)/2 + l;
   double x = ((i - 2 + 0.5)*(1./(1 << point.level))*2. - 0.5);
   double y = ((j - 2 + 0.5)*(1./(1 << point.level))*2. - 0.5);
   if (x > -0.5 && x < 0.5 && y > -0.5 && y < 0.5 &&
       !(aparent(k,l,0).flags & active)) {
     FILE * fp = fopen("check_two_one_loc", "w");
     fprintf (fp,
       "# %d %d\n"
       "%g %g\n%g %g\n",
       k, l,
       (((point.i - 2) + 0.5)*(1./(1 << point.level)) - 0.5),
       (((point.j - 2) + 0.5)*(1./(1 << point.level)) - 0.5),
       x, y);
     fclose (fp);





     if (!(false)) qassert ("/u/basilisk/src/grid/tree.h", 0, "false");
   }
 }end_foreach_leaf();}
}


struct _locate { double x, y, z; };

Point locate (struct _locate p)
{
  for (int l = depth(); l >= 0; l--) {
    Point point = {0};int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
    point.level = l;
    int n = 1 << point.level;
    point.i = (p.x - X0)/L0*n + 2;

    point.j = (p.y - Y0)/L0*n + 2;




    if (point.i >= 0 && point.i < n + 2*2

 && point.j >= 0 && point.j < n + 2*2




 ) {
      if (allocated(0,0,0) && is_local(cell) && is_leaf(cell))
 return point;
    }
    else
      break;
  }
  Point point = {0};int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  point.level = -1;
  return point;
}



bool tree_is_full()
{
  { if (((Tree *)grid)->dirty) update_cache_f(); };
  return (grid->tn == 1L << grid->maxdepth*2);
}

#line 1 "grid/tree-common.h"
#line 1 "/u/basilisk/src/grid/tree-common.h"



#line 1 "grid/multigrid-common.h"
#line 1 "/u/basilisk/src/grid/multigrid-common.h"


#line 1 "grid/cartesian-common.h"
#line 1 "/u/basilisk/src/grid/cartesian-common.h"
#line 1 "grid/events.h"
#line 1 "/u/basilisk/src/grid/events.h"







static void event_error (Event * ev, const char * s)
{
  fprintf (ferr, "%s:%d: error: %s\n", ev->file, ev->line, s);
  exit (1);
}

static void init_event (Event * ev)
{
  if (ev->arrayi || ev->arrayt) {
    ev->i = ev->t = -1;
    if (ev->arrayi)
      ev->i = ev->arrayi[0];
    else
      ev->t = ev->arrayt[0];
    ev->a = 1;
    ev->expr[1] = NULL;
  }
  else {
    if (ev->nexpr > 0) {
      Expr init = NULL, cond = NULL, inc = NULL;
      for (int j = 0; j < ev->nexpr; j++) {
 int i = -123456; double t = i;
 (* ev->expr[j]) (&i, &t, ev);
 if (i == -123456 && t == -123456) {

   if (cond)
     event_error (ev, "events can only use a single condition");
   cond = ev->expr[j];
 }
 else {

   int i1 = i; double t1 = t;
   (* ev->expr[j]) (&i1, &t1, ev);
   if (i1 == i && t1 == t) {


     if (init)
       event_error (ev, "events can only use a single initialisation");
     init = ev->expr[j];
   }
   else {

     if (inc)
       event_error (ev, "events can only use a single increment");
     inc = ev->expr[j];
   }
 }
      }
      ev->expr[0] = init;
      ev->expr[1] = cond;
      ev->expr[2] = inc;
      ev->nexpr = 0;
    }
    ev->i = ev->t = -1;
    if (ev->expr[0]) {
      (* ev->expr[0]) (&ev->i, &ev->t, ev);
      if (ev->i == 1234567890 || ev->t == 1234567890) {
 ev->i = 1234567890; ev->t = -1;
      }
    }
    else if (ev->expr[2]) {
      (* ev->expr[2]) (&ev->i, &ev->t, ev);
      if (ev->i != -1)
 ev->i = 0;
      if (ev->t != -1)
 ev->t = 0;
    }
  }
}

enum { event_done, event_alive, event_stop };

static int event_finished (Event * ev)
{
  ev->t = ev->i = -1;
  return event_done;
}

void event_register (Event event) {
  if (!(Events)) qassert ("/u/basilisk/src/grid/events.h", 0, "Events");
  if (!(!event.last)) qassert ("/u/basilisk/src/grid/events.h", 0, "!event.last");
  int n = 0, parent = -1;
  for (Event * ev = Events; !ev->last; ev++) {
    if (!strcmp (event.name, ev->name)) {
      if (!(parent < 0)) qassert ("/u/basilisk/src/grid/events.h", 0, "parent < 0");
      parent = n;
    }
    n++;
  }
  if (parent < 0) {
    Events = (Event *) prealloc (Events, (n + 2)*sizeof(Event),__func__,__FILE__,0);
    Events[n] = event;
    Events[n].next = NULL;
    Events[n + 1].last = true;
    init_event (&Events[n]);
  }
  else {
    Event * ev = ((Event *) pcalloc (1, sizeof(Event),__func__,__FILE__,0));
    *ev = Events[parent];
    Events[parent] = event;
    Events[parent].next = ev;
    init_event (&Events[parent]);
  }
}

static int event_cond (Event * ev, int i, double t)
{
  if (!ev->expr[1])
    return true;
  return (* ev->expr[1]) (&i, &t, ev);
}
#line 131 "/u/basilisk/src/grid/events.h"
static int event_do (Event * ev, bool action)
{
  if ((iter > ev->i && t > ev->t) || !event_cond (ev, iter, t))
    return event_finished (ev);
  if (iter == ev->i || fabs (t - ev->t) <= 1e-9) {
    if (action) {
      bool finished = false;
      for (Event * e = ev; e; e = e->next) {



 if ((* e->action) (iter, t, e))
   finished = true;
      }
      if (finished) {
 event_finished (ev);
 return event_stop;
      }
    }
    if (ev->arrayi) {
      ev->i = ev->arrayi[ev->a++];
      if (ev->i < 0)
 return event_finished (ev);
    }
    if (ev->arrayt) {
      ev->t = ev->arrayt[ev->a++];
      if (ev->t < 0)
 return event_finished (ev);
    }
    else if (ev->expr[2]) {
      int i0 = ev->i;
      (* ev->expr[2]) (&ev->i, &ev->t, ev);
      if (i0 == -1 && ev->i != i0)
 ev->i += iter + 1;
      if (!event_cond (ev, iter + 1, ev->t))
 return event_finished (ev);
    }
    else if (ev->expr[0] && !ev->expr[1])
      return event_finished (ev);
  }
  return event_alive;
}

static void end_event_do (bool action)
{




  for (Event * ev = Events; !ev->last; ev++)
    if (ev->i == 1234567890 && action)
      for (Event * e = ev; e; e = e->next) {



 e->action (iter, t, e);
      }
}

int events (bool action)
{





  if (iter == 0)
    for (Event * ev = Events; !ev->last; ev++)
      init_event (ev);

  int cond = 0, cond1 = 0;
  inext = 1234567890; tnext = HUGE;
  for (Event * ev = Events; !ev->last && !cond; ev++)
    if (ev->i != 1234567890 &&
 (ev->expr[1] || (ev->expr[0] && !ev->expr[1] && !ev->expr[2]) || ev->arrayi || ev->arrayt))
      cond = 1;
  for (Event * ev = Events; !ev->last; ev++) {
    int status = event_do (ev, action);
    if (status == event_stop) {
      end_event_do (action);
      return 0;
    }
    if (status == event_alive && ev->i != 1234567890 &&
 (ev->expr[1] || (ev->expr[0] && !ev->expr[1] && !ev->expr[2]) || ev->arrayi || ev->arrayt))
      cond1 = 1;
    if (ev->t > t && ev->t < tnext)
      tnext = ev->t;
    if (ev->i > iter && ev->i < inext)
      inext = ev->i;
  }
  if ((!cond || cond1) && (tnext != HUGE || inext != 1234567890)) {
    inext = iter + 1;
    return 1;
  }
  end_event_do (action);
  return 0;
}

void event (const char * name)
{
  for (Event * ev = Events; !ev->last; ev++)
    if (!strcmp (ev->name, name))
      for (Event * e = ev; e; e = e->next) {



 (* e->action) (0, 0, e);
      }
}

double dtnext (double dt)
{
  if (tnext != HUGE && tnext > t) {
    unsigned int n = (tnext - t)/dt;
    if (!(n < INT_MAX)) qassert ("/u/basilisk/src/grid/events.h", 0, "n < INT_MAX");
    if (n == 0)
      dt = tnext - t;
    else {
      double dt1 = (tnext - t)/n;
      if (dt1 > dt + 1e-9)
 dt = (tnext - t)/(n + 1);
      else if (dt1 < dt)
 dt = dt1;
      tnext = t + dt;
    }
  }
  else
    tnext = t + dt;
  return dt;
}
#line 2 "/u/basilisk/src/grid/cartesian-common.h"

void (* debug) (Point);

#define _val_constant(a,k,l,m) ((const double) _constant[a.i -_NVARMAX])
#define diagonalize(a)
#define val_diagonal(a,k,l,m) ((k) == 0 && (l) == 0 && (m) == 0)

#undef VARIABLES
#define VARIABLES\
  double Delta = L0*(1./(1 << point.level));\
  double Delta_x = Delta;\
\
  double Delta_y = Delta;\
\
\
\
\
\
  double x = (ig/2. + (point.i - 2) + 0.5)*Delta + X0; NOT_UNUSED(x);\
\
  double y = (jg/2. + (point.j - 2) + 0.5)*Delta + Y0;\
\
\
\
 NOT_UNUSED(y);\
\
\
\
  double z = 0.;\
\
  NOT_UNUSED(z);\
\
  NOT_UNUSED(Delta);\
  NOT_UNUSED(Delta_x);\
\
  NOT_UNUSED(Delta_y);\
\
\
\
\
\
  ;\

#line 44


#line 1 "grid/fpe.h"
#line 1 "/u/basilisk/src/grid/fpe.h"


#include <signal.h>
#include <unistd.h>

static int gdb()
{
  if (last_point.level >= 0) {
    debug (last_point);
    fputc ('\n', ferr);
    fflush (ferr);
  }
  char command[80];
  sprintf (command, "exec xterm -e 'gdb -p %d' & xterm -e 'gnuplot plot -'",
    getpid());
  return system (command);
}

static void caught_abort (int sig)
{
  fprintf (ferr, "Caught signal %d (Aborted)\n", sig);
  gdb();
}

static void caught_fpe (int sig)
{
  fprintf (ferr, "Caught signal %d (Floating Point Exception)\n", sig);
  gdb();
  exit (1);
}

static void caught_segfault (int sig)
{
  fprintf (ferr, "Caught signal %d (Segmentation Fault)\n", sig);
  gdb();
  exit (2);
}

void catch_fpe (void)
{
  struct sigaction act;
  act.sa_handler = caught_fpe;
  sigemptyset (&act.sa_mask);
  act.sa_flags = 0;
  last_point.level = -1;
  sigaction (8, &act, NULL);
  act.sa_handler = caught_segfault;
  sigaction (11, &act, NULL);
  act.sa_handler = caught_abort;
  act.sa_flags = SA_RESETHAND;
  sigaction (6, &act, NULL);
}
#line 47 "/u/basilisk/src/grid/cartesian-common.h"

#define end_foreach_face()



static void init_block_scalar (scalar sb, const char * name, const char * ext,
          int n, int block)
{
  char bname[strlen(name) + strlen(ext) + 10];
  if (n == 0) {
    sprintf (bname, "%s%s", name, ext);
    _attribute[sb.i].block = block;
    init_scalar (sb, bname);
    baseblock = list_append (baseblock, sb);
  }
  else {
    sprintf (bname, "%s%d%s", name, n, ext);
    _attribute[sb.i].block = - n;
    init_scalar (sb, bname);
  }
  all = list_append (all, sb);
}

scalar new_block_scalar (const char * name, const char * ext, int block)
{
  int nvar = datasize/sizeof(double);

  scalar s = {0};
  while (s.i < nvar) {
    int n = 0;
    scalar sb = s;
    while (sb.i < nvar && n < block && _attribute[sb.i].freed)
      n++, sb.i++;
    if (n >= block) {
      for (sb.i = s.i, n = 0; n < block; n++, sb.i++)
 init_block_scalar (sb, name, ext, n, block);
      trash (((scalar []){s, {-1}}));
      return s;
    }
    s.i = sb.i + 1;
  }


  s = (scalar){nvar};
  if (!(nvar + block <= _NVARMAX)) qassert ("/u/basilisk/src/grid/cartesian-common.h", 0, "nvar + block <= _NVARMAX");
  _attribute = (_Attributes *) prealloc (_attribute, (nvar + block)*sizeof(_Attributes),__func__,__FILE__,0);
  memset (&_attribute[nvar], 0, block*sizeof (_Attributes));
  for (int n = 0; n < block; n++, nvar++) {
    scalar sb = (scalar){nvar};
    init_block_scalar (sb, name, ext, n, block);
  }

  realloc_scalar (block*sizeof(double));
  trash (((scalar []){s, {-1}}));
  return s;
}

scalar new_scalar (const char * name)
{
  return new_block_scalar (name, "", 1);
}

scalar new_vertex_scalar (const char * name)
{
  scalar s = new_block_scalar (name, "", 1);
  init_vertex_scalar (s, NULL);
  return s;
}

scalar new_block_vertex_scalar (const char * name, int block)
{
  scalar s = new_block_scalar (name, "", block);
  for (int i = 0; i < block; i++) {
    scalar sb = {s.i + i};
    init_vertex_scalar (sb, NULL);
  }
  return s;
}

static vector alloc_block_vector (const char * name, int block)
{
  vector v;
  struct { char * x, * y, * z; } ext = {".x", ".y", ".z"};
  
    v.x = new_block_scalar (name, ext.x, block);
    
#line 131
v.y = new_block_scalar (name, ext.y, block);
  return v;
}

vector new_vector (const char * name)
{
  vector v = alloc_block_vector (name, 1);
  init_vector (v, NULL);
  return v;
}

vector new_face_vector (const char * name)
{
  vector v = alloc_block_vector (name, 1);
  init_face_vector (v, NULL);
  return v;
}

vector new_block_vector (const char * name, int block)
{
  vector v = alloc_block_vector (name, block);
  for (int i = 0; i < block; i++) {
    vector vb;
    
      vb.x.i = v.x.i + i;
      
#line 155
vb.y.i = v.y.i + i;
    init_vector (vb, NULL);
    
      _attribute[vb.x.i].block = - i;
      
#line 158
_attribute[vb.y.i].block = - i;
  }
  
    _attribute[v.x.i].block = block;
    
#line 161
_attribute[v.y.i].block = block;
  return v;
}

vector new_block_face_vector (const char * name, int block)
{
  vector v = alloc_block_vector (name, block);
  for (int i = 0; i < block; i++) {
    vector vb;
    
      vb.x.i = v.x.i + i;
      
#line 171
vb.y.i = v.y.i + i;
    init_face_vector (vb, NULL);
    
      _attribute[vb.x.i].block = - i;
      
#line 174
_attribute[vb.y.i].block = - i;
  }
  
    _attribute[v.x.i].block = block;
    
#line 177
_attribute[v.y.i].block = block;
  return v;
}

tensor new_tensor (const char * name)
{
  char cname[strlen(name) + 3];
  struct { char * x, * y, * z; } ext = {"%s.x", "%s.y", "%s.z"};
  tensor t;
   {
    sprintf (cname, ext.x, name);
    t.x = new_vector (cname);
  } 
#line 186
{
    sprintf (cname, ext.y, name);
    t.y = new_vector (cname);
  }
  init_tensor (t, NULL);
  return t;
}

tensor new_symmetric_tensor (const char * name)
{
  char cname[strlen(name) + 5];
  struct { char * x, * y, * z; } ext = {"%s.x.x", "%s.y.y", "%s.z.z"};
  tensor t;
   {
    sprintf (cname, ext.x, name);
    t.x.x = new_scalar(cname);
  } 
#line 199
{
    sprintf (cname, ext.y, name);
    t.y.y = new_scalar(cname);
  }

    sprintf (cname, "%s.x.y", name);
    t.x.y = new_scalar(cname);
    t.y.x = t.x.y;
#line 219 "/u/basilisk/src/grid/cartesian-common.h"
  init_tensor (t, NULL);
  return t;
}

static int nconst = 0;

void init_const_scalar (scalar s, const char * name, double val)
{
  if (s.i - _NVARMAX >= nconst) {
    nconst = s.i - _NVARMAX + 1;
    _constant = (double *) prealloc (_constant, (nconst)*sizeof(double),__func__,__FILE__,0);
  }
  _constant[s.i - _NVARMAX] = val;
}

scalar new_const_scalar (const char * name, int i, double val)
{
  scalar s = (scalar){i + _NVARMAX};
  init_const_scalar (s, name, val);
  return s;
}

void init_const_vector (vector v, const char * name, double * val)
{
  
    init_const_scalar (v.x, name, *val++);
    
#line 244
init_const_scalar (v.y, name, *val++);
}

vector new_const_vector (const char * name, int i, double * val)
{
  vector v;
  
    v.x.i = _NVARMAX + i++;
    
#line 251
v.y.i = _NVARMAX + i++;
  init_const_vector (v, name, val);
  return v;
}

void scalar_clone (scalar a, scalar b)
{
  char * name = _attribute[a.i].name;
  double (** boundary) (Point, Point, scalar, void *) = _attribute[a.i].boundary;
  double (** boundary_homogeneous) (Point, Point, scalar, void *) =
    _attribute[a.i].boundary_homogeneous;
  if (!(_attribute[b.i].block > 0 && _attribute[a.i].block == _attribute[b.i].block)) qassert ("/u/basilisk/src/grid/cartesian-common.h", 0, "b.block > 0 && a.block == b.block");
  pfree (_attribute[a.i].depends,__func__,__FILE__,0);
  _attribute[a.i] = _attribute[b.i];
  _attribute[a.i].name = name;
  _attribute[a.i].boundary = boundary;
  _attribute[a.i].boundary_homogeneous = boundary_homogeneous;
  for (int i = 0; i < nboundary; i++) {
    _attribute[a.i].boundary[i] = _attribute[b.i].boundary[i];
    _attribute[a.i].boundary_homogeneous[i] = _attribute[b.i].boundary_homogeneous[i];
  }
  _attribute[a.i].depends = list_copy (_attribute[b.i].depends);
}

scalar * list_clone (scalar * l)
{
  scalar * list = NULL;
  int nvar = datasize/sizeof(double), map[nvar];
  for (int i = 0; i < nvar; i++)
    map[i] = -1;
  {scalar*_i=(scalar*)( l);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {
    scalar c = _attribute[s.i].block > 1 ? new_block_scalar("c", "", _attribute[s.i].block) :
      new_scalar("c");
    scalar_clone (c, s);
    map[s.i] = c.i;
    list = list_append (list, c);
  }}}
  {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    {
      if (_attribute[s.i].v.x.i >= 0 && map[_attribute[s.i].v.x.i] >= 0)
 _attribute[s.i].v.x.i = map[_attribute[s.i].v.x.i];
      
#line 290
if (_attribute[s.i].v.y.i >= 0 && map[_attribute[s.i].v.y.i] >= 0)
 _attribute[s.i].v.y.i = map[_attribute[s.i].v.y.i];}}}
  return list;
}

void delete (scalar * list)
{
  if (all == NULL)
    return;

  {scalar*_i=(scalar*)( list);if(_i)for(scalar f=*_i;(&f)->i>=0;f=*++_i){ {
    for (int i = 0; i < _attribute[f.i].block; i++) {
      scalar fb = {f.i + i};
      if (_attribute[f.i].delete)
 _attribute[f.i].delete (fb);
      pfree (_attribute[fb.i].name,__func__,__FILE__,0); _attribute[fb.i].name = NULL;
      pfree (_attribute[fb.i].boundary,__func__,__FILE__,0); _attribute[fb.i].boundary = NULL;
      pfree (_attribute[fb.i].boundary_homogeneous,__func__,__FILE__,0); _attribute[fb.i].boundary_homogeneous = NULL;
      pfree (_attribute[fb.i].depends,__func__,__FILE__,0); _attribute[fb.i].depends = NULL;
      _attribute[fb.i].freed = true;
    }
  }}}

  if (list == all) {
    all[0].i = -1;
    baseblock[0].i = -1;
    return;
  }

  trash (list);
  {scalar*_i=(scalar*)( list);if(_i)for(scalar f=*_i;(&f)->i>=0;f=*++_i){ {
    if (_attribute[f.i].block > 0) {
      scalar * s = all;
      for (; s->i >= 0 && s->i != f.i; s++);
      if (s->i == f.i) {
 for (; s[_attribute[f.i].block].i >= 0; s++)
   s[0] = s[_attribute[f.i].block];
 s->i = -1;
      }
      for (s = baseblock; s->i >= 0 && s->i != f.i; s++);
      if (s->i == f.i) {
 for (; s[1].i >= 0; s++)
   s[0] = s[1];
 s->i = -1;
      }
    }
  }}}
}

void free_solver()
{
  if (!(_val_higher_dimension == 0.)) qassert ("/u/basilisk/src/grid/cartesian-common.h", 0, "_val_higher_dimension == 0.");

  if (free_solver_funcs) {
    free_solver_func * a = (free_solver_func *) free_solver_funcs->p;
    for (int i = 0; i < free_solver_funcs->len/sizeof(free_solver_func); i++)
      a[i] ();
    array_free (free_solver_funcs);
  }

  delete (all);
  pfree (all,__func__,__FILE__,0); all = NULL;
  pfree (baseblock,__func__,__FILE__,0); baseblock = NULL;
  for (Event * ev = Events; !ev->last; ev++) {
    Event * e = ev->next;
    while (e) {
      Event * next = e->next;
      pfree (e,__func__,__FILE__,0);
      e = next;
    }
  }

  pfree (Events,__func__,__FILE__,0); Events = NULL;
  pfree (_attribute,__func__,__FILE__,0); _attribute = NULL;
  pfree (_constant,__func__,__FILE__,0); _constant = NULL;
  free_grid();
  qpclose_all();
#if TRACE
  trace_off();
#endif
#if MTRACE
  pmuntrace();
#endif
#if _CADNA
  cadna_end();
#endif
}



void (* boundary_level) (scalar *, int l);
void (* boundary_face) (vectorl);




void boundary_flux (vector * list) __attribute__ ((deprecated));

void boundary_flux (vector * list)
{
  vectorl list1 = {NULL};
  {vector*_i=(vector*)( list);if(_i)for(vector v=*_i;(&v)->x.i>=0;v=*++_i){
    {
      list1.x = list_append (list1.x, v.x);
      
#line 393
list1.y = list_append (list1.y, v.y);}}}
  boundary_face (list1);
  
    pfree (list1.x,__func__,__FILE__,0);
    
#line 396
pfree (list1.y,__func__,__FILE__,0);
}

static scalar * list_add_depends (scalar * list, scalar s)
{
  {scalar*_i=(scalar*)( list);if(_i)for(scalar t=*_i;(&t)->i>=0;t=*++_i){
    if (t.i == s.i)
      return list;}}
  scalar * list1 = list;
  {scalar*_i=(scalar*)( _attribute[s.i].depends);if(_i)for(scalar d=*_i;(&d)->i>=0;d=*++_i){
    if (_attribute[d.i].dirty)
      list1 = list_add_depends (list1, d);}}
  return list_append (list1, s);
}

     
void boundary_internal (scalar * list, const char * fname, int line)
{tracing("boundary_internal","/u/basilisk/src/grid/cartesian-common.h",0);
  if (list == NULL)
    {end_tracing("boundary_internal","/u/basilisk/src/grid/cartesian-common.h",0);return;}
  scalar * listc = NULL;
  vectorl listf = {NULL};
  bool flux = false;
  {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    if (!is_constant(s) && _attribute[s.i].block > 0) {
      if (scalar_is_dirty (s)) {
 if (_attribute[s.i].face && _attribute[s.i].dirty != 2)
   {
     if (_attribute[s.i].v.x.i == s.i)
       listf.x = list_add (listf.x, s), flux = true;
     
#line 424
if (_attribute[s.i].v.y.i == s.i)
       listf.y = list_add (listf.y, s), flux = true;}
 if (!is_constant(cm) && _attribute[cm.i].dirty)
   listc = list_add_depends (listc, cm);
 if (_attribute[s.i].face != 2)
   listc = list_add_depends (listc, s);
      }




    }}}
  if (flux) {
    boundary_face (listf);
    
      pfree (listf.x,__func__,__FILE__,0);
      
#line 439
pfree (listf.y,__func__,__FILE__,0);
  }
  if (listc) {
    boundary_level (listc, -1);
    {scalar*_i=(scalar*)( listc);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      _attribute[s.i].dirty = false;}}
    pfree (listc,__func__,__FILE__,0);
  }
end_tracing("boundary_internal","/u/basilisk/src/grid/cartesian-common.h",0);}

void cartesian_boundary_level (scalar * list, int l)
{
  { Boundary ** _i = boundaries, * _b; while (_i && (_b = *_i++)) if (_b->level) _b->level (_b, list, l); };
}

void cartesian_boundary_face (vectorl list)
{
  
    {scalar*_i=(scalar*)( list.x);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      _attribute[s.i].dirty = 2;}}
    
#line 457
{scalar*_i=(scalar*)( list.y);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      _attribute[s.i].dirty = 2;}}
}

static double symmetry (Point point, Point neighbor, scalar s, void * data)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  return val(s,0,0,0);
}

static double antisymmetry (Point point, Point neighbor, scalar s, void * data)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  return -val(s,0,0,0);
}

double (* default_scalar_bc[]) (Point, Point, scalar, void *) = {
  symmetry, symmetry, symmetry, symmetry, symmetry, symmetry
};

scalar cartesian_init_scalar (scalar s, const char * name)
{

  char * pname;
  if (name) {
    pfree (_attribute[s.i].name,__func__,__FILE__,0);
    pname = pstrdup (name,__func__,__FILE__,0);
  }
  else
    pname = _attribute[s.i].name;
  int block = _attribute[s.i].block;
  double (** boundary) (Point, Point, scalar, void *) = _attribute[s.i].boundary;
  double (** boundary_homogeneous) (Point, Point, scalar, void *) =
    _attribute[s.i].boundary_homogeneous;

  _attribute[s.i] = (const _Attributes){0};
  _attribute[s.i].name = pname;
  _attribute[s.i].block = block == 0 ? 1 : block;

  _attribute[s.i].boundary = boundary ? boundary :
    (double (**)(Point, Point, scalar, void *))
    pmalloc (nboundary*sizeof (void (*)()),__func__,__FILE__,0);
  _attribute[s.i].boundary_homogeneous = boundary_homogeneous ? boundary_homogeneous :
    (double (**)(Point, Point, scalar, void *))
    pmalloc (nboundary*sizeof (void (*)()),__func__,__FILE__,0);
  for (int b = 0; b < nboundary; b++)
    _attribute[s.i].boundary[b] = _attribute[s.i].boundary_homogeneous[b] =
      b < 2*2 ? default_scalar_bc[b] : symmetry;
  _attribute[s.i].gradient = NULL;
   {
    _attribute[s.i].d.x = 0;
    _attribute[s.i].v.x.i = -1;
  } 
#line 504
{
    _attribute[s.i].d.y = 0;
    _attribute[s.i].v.y.i = -1;
  }
  _attribute[s.i].face = false;
  return s;
}

scalar cartesian_init_vertex_scalar (scalar s, const char * name)
{
  s = cartesian_init_scalar (s, name);
  
    _attribute[s.i].d.x = -1;
    
#line 516
_attribute[s.i].d.y = -1;
  for (int d = 0; d < nboundary; d++)
    _attribute[s.i].boundary[d] = _attribute[s.i].boundary_homogeneous[d] = NULL;
  return s;
}

double (* default_vector_bc[]) (Point, Point, scalar, void *) = {
  antisymmetry, antisymmetry,
  antisymmetry, antisymmetry,
  antisymmetry, antisymmetry
};

vector cartesian_init_vector (vector v, const char * name)
{
  struct { char * x, * y, * z; } ext = {".x", ".y", ".z"};
   {
    if (name) {
      char cname[strlen(name) + 3];
      sprintf (cname, "%s%s", name, ext.x);
      init_scalar (v.x, cname);
    }
    else
      init_scalar (v.x, NULL);
    _attribute[v.x.i].v = v;
  } 
#line 531
{
    if (name) {
      char cname[strlen(name) + 3];
      sprintf (cname, "%s%s", name, ext.y);
      init_scalar (v.y, cname);
    }
    else
      init_scalar (v.y, NULL);
    _attribute[v.y.i].v = v;
  }

  for (int d = 0; d < nboundary; d++)
    _attribute[v.x.i].boundary[d] = _attribute[v.x.i].boundary_homogeneous[d] =
      d < 2*2 ? default_vector_bc[d] : antisymmetry;
  return v;
}

vector cartesian_init_face_vector (vector v, const char * name)
{
  v = cartesian_init_vector (v, name);
   {
    _attribute[v.x.i].d.x = -1;
    _attribute[v.x.i].face = true;
  } 
#line 551
{
    _attribute[v.y.i].d.y = -1;
    _attribute[v.y.i].face = true;
  }
  for (int d = 0; d < nboundary; d++)
    _attribute[v.x.i].boundary[d] = _attribute[v.x.i].boundary_homogeneous[d] = NULL;
  return v;
}

tensor cartesian_init_tensor (tensor t, const char * name)
{
  struct { char * x, * y, * z; } ext = {".x", ".y", ".z"};
   {
    if (name) {
      char cname[strlen(name) + 3];
      sprintf (cname, "%s%s", name, ext.x);
      init_vector (t.x, cname);
    }
    else
      init_vector (t.x, NULL);
  } 
#line 563
{
    if (name) {
      char cname[strlen(name) + 3];
      sprintf (cname, "%s%s", name, ext.y);
      init_vector (t.y, cname);
    }
    else
      init_vector (t.y, NULL);
  }






    for (int b = 0; b < nboundary; b++) {
      _attribute[t.x.x.i].boundary[b] = _attribute[t.y.x.i].boundary[b] =
 _attribute[t.x.x.i].boundary_homogeneous[b] = _attribute[t.y.y.i].boundary_homogeneous[b] =
 b < 2*2 ? default_scalar_bc[b] : symmetry;
      _attribute[t.x.y.i].boundary[b] = _attribute[t.y.y.i].boundary[b] =
 _attribute[t.x.y.i].boundary_homogeneous[b] = _attribute[t.y.x.i].boundary_homogeneous[b] =
 b < 2*2 ? default_vector_bc[b] : antisymmetry;
    }



  return t;
}

struct OutputCells {
  FILE * fp;
  coord c;
  double size;
};

void output_cells (struct OutputCells p)
{
  if (!p.fp) p.fp = fout;
  {foreach() {
    bool inside = true;
    coord o = {x,y,z};
    
      if (inside && p.size > 0. &&
   (o.x > p.c.x + p.size || o.x < p.c.x - p.size))
 inside = false;
      
#line 605
if (inside && p.size > 0. &&
   (o.y > p.c.y + p.size || o.y < p.c.y - p.size))
 inside = false;
    if (inside) {
      Delta /= 2.;



      fprintf (p.fp, "%g %g\n%g %g\n%g %g\n%g %g\n%g %g\n\n",
        x - Delta, y - Delta,
        x - Delta, y + Delta,
        x + Delta, y + Delta,
        x + Delta, y - Delta,
        x - Delta, y - Delta);
#line 633 "/u/basilisk/src/grid/cartesian-common.h"
    }
  }end_foreach();}
  fflush (p.fp);
}
#line 645 "/u/basilisk/src/grid/cartesian-common.h"
static char * replace_ (const char * vname)
{
  char * name = pstrdup (vname,__func__,__FILE__,0), * c = name;
  while (*c != '\0') {
    if (*c == '.')
      *c = '_';
    c++;
  }
  return name;
}

static void debug_plot (FILE * fp, const char * name, const char * cells,
   const char * stencil)
{
  char * vname = replace_ (name);
  fprintf (fp,
    "  load 'debug.plot'\n"
    "  v=%s\n"




    "  plot '%s' w l lc 0, "
    "'%s' u 1+3*v:2+3*v:3+3*v w labels tc lt 1 title columnhead(3+3*v)",





    vname, cells, stencil);
  pfree (vname,__func__,__FILE__,0);
}

void cartesian_debug (Point point)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  char name[80] = "cells";
  if (pid() > 0)
    sprintf (name, "cells-%d", pid());
  FILE * fp = fopen (name, "w");
  output_cells ((struct OutputCells){fp, (coord){x,y,z}, 4.*Delta});
  fclose (fp);

  char stencil[80] = "stencil";
  if (pid() > 0)
    sprintf (stencil, "stencil-%d", pid());
  fp = fopen (stencil, "w");
  {scalar*_i=(scalar*)( all);if(_i)for(scalar v=*_i;(&v)->i>=0;v=*++_i){



    fprintf (fp, "x y %s ", _attribute[v.i].name);}}



  fputc ('\n', fp);
#line 712 "/u/basilisk/src/grid/cartesian-common.h"
    for (int k = -2; k <= 2; k++)
      for (int l = -2; l <= 2; l++) {
 {scalar*_i=(scalar*)( all);if(_i)for(scalar v=*_i;(&v)->i>=0;v=*++_i){ {
   fprintf (fp, "%g %g ",
     x + k*Delta + _attribute[v.i].d.x*Delta/2.,
     y + l*Delta + _attribute[v.i].d.y*Delta/2.);
   if (allocated(k,l,0))
     fprintf (fp, "%g ", val(v,k,l,0));
   else
     fputs ("n/a ", fp);
 }}}
 fputc ('\n', fp);
      }
#line 742 "/u/basilisk/src/grid/cartesian-common.h"
  fclose (fp);

  fp = fopen ("debug.plot", "w");
  fprintf (fp,
    "set term x11\n"
    "set size ratio -1\n"
    "set key outside\n");
  {scalar*_i=(scalar*)( all);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {
    char * name = replace_ (_attribute[s.i].name);
    fprintf (fp, "%s = %d\n", name, s.i);
    pfree (name,__func__,__FILE__,0);
  }}}
  fclose (fp);

  fprintf (ferr, "Last point stencils can be displayed using (in gnuplot)\n");
  debug_plot (ferr, _attribute[0].name, name, stencil);
  fflush (ferr);

  fp = fopen ("plot", "w");
  debug_plot (fp, _attribute[0].name, name, stencil);
  fclose (fp);
}

void cartesian_methods()
{
  init_scalar = cartesian_init_scalar;
  init_vertex_scalar = cartesian_init_vertex_scalar;
  init_vector = cartesian_init_vector;
  init_tensor = cartesian_init_tensor;
  init_face_vector = cartesian_init_face_vector;
  boundary_level = cartesian_boundary_level;
  boundary_face = cartesian_boundary_face;
  debug = cartesian_debug;
}

tensor init_symmetric_tensor (tensor t, const char * name)
{
  return init_tensor (t, name);
}

struct _interpolate {
  scalar v;
  double x, y, z;
};

static double interpolate_linear (Point point, struct _interpolate p)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  scalar v = p.v;







  x = (p.x - x)/Delta - _attribute[v.i].d.x/2.;
  y = (p.y - y)/Delta - _attribute[v.i].d.y/2.;
  int i = sign(x), j = sign(y);
  x = fabs(x); y = fabs(y);

  return ((val(v,0,0,0)*(1. - x) + val(v,i,0,0)*x)*(1. - y) +
   (val(v,0,j,0)*(1. - x) + val(v,i,j,0)*x)*y);
#line 816 "/u/basilisk/src/grid/cartesian-common.h"
}

     
double interpolate (struct _interpolate p)
{tracing("interpolate","/u/basilisk/src/grid/cartesian-common.h",0);
  scalar v = p.v;
  boundary_internal ((scalar *)((scalar[]){v,{-1}}), "/u/basilisk/src/grid/cartesian-common.h", 0);
  Point point = locate ((struct _locate){p.x, p.y, p.z});int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  if (point.level < 0)
    {end_tracing("interpolate","/u/basilisk/src/grid/cartesian-common.h",0);return HUGE;}
  { double _ret= interpolate_linear (point, p);end_tracing("interpolate","/u/basilisk/src/grid/cartesian-common.h",0);return _ret;}
end_tracing("interpolate","/u/basilisk/src/grid/cartesian-common.h",0);}

     
void interpolate_array (scalar * list, coord * a, int n, double * v, bool linear)
{tracing("interpolate_array","/u/basilisk/src/grid/cartesian-common.h",0);
  boundary_internal ((scalar *)list, "/u/basilisk/src/grid/cartesian-common.h", 0);
  int j = 0;
  for (int i = 0; i < n; i++) {
    Point point = locate ((struct _locate){a[i].x, a[i].y, a[i].z});int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
    if (point.level >= 0) {
      {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
 v[j++] = !linear ? val(s,0,0,0) :
   interpolate_linear (point,
         (struct _interpolate){s, a[i].x, a[i].y, a[i].z});}}
    }
    else
      {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
 v[j++] = HUGE;}}
  }
#if _MPI
  if (pid() == 0)
    MPI_Reduce (MPI_IN_PLACE, v, n*list_len(list), MPI_DOUBLE,
  MPI_MIN, 0, MPI_COMM_WORLD);
  else
    MPI_Reduce (v, v, n*list_len(list), MPI_DOUBLE,
  MPI_MIN, 0, MPI_COMM_WORLD);
#endif
end_tracing("interpolate_array","/u/basilisk/src/grid/cartesian-common.h",0);}



typedef int bid;

bid new_bid()
{
  int b = nboundary++;
  {scalar*_i=(scalar*)( all);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {
    _attribute[s.i].boundary = (double (**)(Point, Point, scalar, void *))
      prealloc (_attribute[s.i].boundary, nboundary*sizeof (void (*)()),__func__,__FILE__,0);
    _attribute[s.i].boundary_homogeneous = (double (**)(Point, Point, scalar, void *))
      prealloc (_attribute[s.i].boundary_homogeneous, nboundary*sizeof (void (*)()),__func__,__FILE__,0);
  }}}
  {scalar*_i=(scalar*)( all);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {
    if (_attribute[s.i].v.x.i < 0)
      _attribute[s.i].boundary[b] = _attribute[s.i].boundary_homogeneous[b] = symmetry;
    else if (_attribute[s.i].v.x.i == s.i) {
      vector v = _attribute[s.i].v;
      
 _attribute[v.y.i].boundary[b] = _attribute[v.y.i].boundary_homogeneous[b] = symmetry;
 
#line 875
_attribute[v.x.i].boundary[b] = _attribute[v.x.i].boundary_homogeneous[b] = symmetry;
      _attribute[v.x.i].boundary[b] = _attribute[v.x.i].boundary_homogeneous[b] =
 _attribute[v.x.i].face ? NULL : antisymmetry;
    }
  }}}
  return b;
}



static double periodic_bc (Point point, Point neighbor, scalar s, void * data)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  return HUGE;
}

static void periodic_boundary (int d)
{

  {scalar*_i=(scalar*)( all);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    _attribute[s.i].boundary[d] = _attribute[s.i].boundary_homogeneous[d] = periodic_bc;}}

  {scalar*_i=(scalar*)( all);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    if (_attribute[s.i].face) {
      vector v = _attribute[s.i].v;
      _attribute[v.x.i].boundary[d] = _attribute[v.x.i].boundary_homogeneous[d] = NULL;
    }}}

  default_scalar_bc[d] = periodic_bc;
  default_vector_bc[d] = periodic_bc;
}

void periodic (int dir)
{



    if (!(dir <= bottom)) qassert ("/u/basilisk/src/grid/cartesian-common.h", 0, "dir <= bottom");




  int c = dir/2;
  periodic_boundary (2*c);
  periodic_boundary (2*c + 1);
  (&Period.x)[c] = true;
}


double getvalue (Point point, scalar s, int i, int j, int k)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  return val(s,i,j,k);
}

void default_stencil (Point p, scalar * list)
{
  {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    _attribute[s.i].input = true, _attribute[s.i].width = 2;}}
}




static void write_stencil_index (int * index)
{
  fprintf (qstderr(), "[%d", index[0]);
  for (int d = 1; d < 2; d++)
    fprintf (qstderr(), ",%d", index[d]);
  fputs ("]", qstderr());
}

void stencil_val (Point p, scalar s, int i, int j, int k,
    const char * file, int line, bool overflow)
{
  if (is_constant(s) || s.i < 0)
    return;
  int index[] = {i, j, k};
  for (int d = 0; d < 2; d++)
    index[d] += (&p.i)[d];
  bool central = true;
  for (int d = 0; d < 2; d++) {
    if (!overflow && (index[d] > 2 || index[d] < - 2)) {
      fprintf (qstderr(), "%s:%d: error: stencil overflow: %s",
        file, line, _attribute[s.i].name);
      write_stencil_index (index);
      fprintf (qstderr(), "\n");
      fflush (qstderr());
      abort();
    }
    if (index[d] != 0)
      central = false;
  }
  if (central) {
    if (!_attribute[s.i].output)
      _attribute[s.i].input = true;
  }
  else {
    _attribute[s.i].input = true;
    int d = 0;
     {
      if ((!_attribute[s.i].face || _attribute[s.i].v.x.i != s.i) && abs(index[d]) > _attribute[s.i].width)
 _attribute[s.i].width = abs(index[d]);
      d++;
    } 
#line 973
{
      if ((!_attribute[s.i].face || _attribute[s.i].v.y.i != s.i) && abs(index[d]) > _attribute[s.i].width)
 _attribute[s.i].width = abs(index[d]);
      d++;
    }
  }
}

void stencil_val_a (Point p, scalar s, int i, int j, int k, bool input,
      const char * file, int line)
{
  if (is_constant(s) || s.i < 0)
    abort();
  int index[] = {i, j, k};
  for (int d = 0; d < 2; d++)
    index[d] += (&p.i)[d];
  for (int d = 0; d < 2; d++)
    if (index[d] != 0) {
      fprintf (qstderr(), "%s:%d: error: illegal write: %s",
        file, line, _attribute[s.i].name);
      write_stencil_index (index);
      fprintf (qstderr(), "\n");
      fflush (qstderr());
      abort();
    }
  if (input && !_attribute[s.i].output)
    _attribute[s.i].input = true;
  _attribute[s.i].output = true;
}
#line 4 "/u/basilisk/src/grid/multigrid-common.h"

#ifndef foreach_level_or_leaf
# define foreach_level_or_leaf foreach_level
# define end_foreach_level_or_leaf end_foreach_level
#endif

#ifndef foreach_coarse_level
# define foreach_coarse_level foreach_level
# define end_foreach_coarse_level end_foreach_level
#endif










void (* restriction) (scalar *);

static inline void restriction_average (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  double sum = 0.;
  {foreach_child()
    sum += val(s,0,0,0);end_foreach_child()}
  val(s,0,0,0) = sum/(1 << 2);
}

static inline void restriction_volume_average (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;

#line 35
if(!is_constant(cm)){{
  double sum = 0.;
  {foreach_child()
    sum += val(cm,0,0,0)*val(s,0,0,0);end_foreach_child()}
  val(s,0,0,0) = sum/(1 << 2)/(val(cm,0,0,0) + 1e-30);
}}else {double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);

#line 35
{
  double sum = 0.;
  {foreach_child()
    sum += _const_cm*val(s,0,0,0);end_foreach_child()}
  val(s,0,0,0) = sum/(1 << 2)/(_const_cm + 1e-30);
}}

#line 40
}

static inline void face_average (Point point, vector v)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
   {




      val(v.x,0,0,0) = (fine(v.x,0,0,0) + fine(v.x,0,1,0))/2.;
      val(v.x,1,0,0) = (fine(v.x,2,0,0) + fine(v.x,2,1,0))/2.;






  } 
#line 44
{




      val(v.y,0,0,0) = (fine(v.y,0,0,0) + fine(v.y,1,0,0))/2.;
      val(v.y,0,1,0) = (fine(v.y,0,2,0) + fine(v.y,1,2,0))/2.;






  }
}

static inline void restriction_face (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  face_average (point, _attribute[s.i].v);
}

static inline void restriction_vertex (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  for (int i = 0; i <= 1; i++) {
    val(s,i,0,0) = fine(s,2*i,0,0);

    val(s,i,1,0) = fine(s,2*i,2,0);





  }
}

static inline void no_restriction (Point point, scalar s) {int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;}

static inline void no_data (Point point, scalar s) {int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  {foreach_child()
    val(s,0,0,0) = HUGE;end_foreach_child()}
}

void wavelet (scalar s, scalar w)
{
  restriction (((scalar[]){s,{-1}}));
  for (int l = grid->maxdepth - 1; l >= 0; l--) {
    {foreach_coarse_level (l) {
      {foreach_child()
        val(w,0,0,0) = val(s,0,0,0);end_foreach_child()}
      _attribute[s.i].prolongation (point, s);
      {foreach_child() {
        double sp = val(s,0,0,0);
        val(s,0,0,0) = val(w,0,0,0);

        val(w,0,0,0) -= sp;
      }end_foreach_child()}
    }end_foreach_coarse_level();}
    boundary_level (((scalar[]){w,{-1}}), l + 1);
  }

  {foreach_level(0)
    val(w,0,0,0) = val(s,0,0,0);end_foreach_level();}
  boundary_level (((scalar[]){w,{-1}}), 0);
}

void inverse_wavelet (scalar s, scalar w)
{
  {foreach_level(0)
    val(s,0,0,0) = val(w,0,0,0);end_foreach_level();}
  boundary_level (((scalar[]){s,{-1}}), 0);
  for (int l = 0; l <= grid->maxdepth - 1; l++) {
    {foreach_coarse_level (l) {
      _attribute[s.i].prolongation (point, s);
      {foreach_child()
        val(s,0,0,0) += val(w,0,0,0);end_foreach_child()}
    }end_foreach_coarse_level();}
    boundary_level (((scalar[]){s,{-1}}), l + 1);
  }
}

static inline double bilinear (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;



    return (9.*coarse(s,0,0,0) +
     3.*(coarse(s,child.x,0,0) + coarse(s,0,child.y,0)) +
     coarse(s,child.x,child.y,0))/16.;
#line 140 "/u/basilisk/src/grid/multigrid-common.h"
}

static inline void refine_bilinear (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  {foreach_child()
    val(s,0,0,0) = bilinear (point, s);end_foreach_child()}
}

static inline double quadratic (double a, double b, double c)
{
  return (30.*a + 5.*b - 3.*c)/32.;
}

static inline double biquadratic (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;



  return
    quadratic (quadratic (coarse(s,0,0,0),
     coarse(s,child.x,0,0),
     coarse(s,-child.x,0,0)),
        quadratic (coarse(s,0,child.y,0),
     coarse(s,child.x,child.y,0),
     coarse(s,-child.x,child.y,0)),
        quadratic (coarse(s,0,-child.y,0),
     coarse(s,child.x,-child.y,0),
     coarse(s,-child.x,-child.y,0)));




}

static inline double biquadratic_vertex (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;



  return (36.*val(s,0,0,0) + 18.*(val(s,-1,0,0) + val(s,0,-1,0)) - 6.*(val(s,1,0,0) + val(s,0,1,0)) +
   9.*val(s,-1,-1,0) - 3.*(val(s,1,-1,0) + val(s,-1,1,0)) + val(s,1,1,0))/64.;




}

static inline void refine_biquadratic (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  {foreach_child()
    val(s,0,0,0) = biquadratic (point, s);end_foreach_child()}
}

static inline void refine_linear (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;

#line 194
if(!is_constant(cm)){{
  coord g;
  if (_attribute[s.i].gradient)
    {
      g.x = _attribute[s.i].gradient (val(s,-1,0,0), val(s,0,0,0), val(s,1,0,0));
      
#line 198
g.y = _attribute[s.i].gradient (val(s,0,-1,0), val(s,0,0,0), val(s,0,1,0));}
  else
    {
      g.x = (val(s,1,0,0) - val(s,-1,0,0))/2.;
      
#line 201
g.y = (val(s,0,1,0) - val(s,0,-1,0))/2.;}

  double sc = val(s,0,0,0), cmc = 4.*val(cm,0,0,0), sum = val(cm,0,0,0)*(1 << 2);
  {foreach_child() {
    val(s,0,0,0) = sc;
    
      val(s,0,0,0) += child.x*g.x*val(cm,-child.x,0,0)/cmc;
      
#line 207
val(s,0,0,0) += child.y*g.y*val(cm,0,-child.y,0)/cmc;
    sum -= val(cm,0,0,0);
  }end_foreach_child()}
  if (!(fabs(sum) < 1e-10)) qassert ("/u/basilisk/src/grid/multigrid-common.h", 0, "fabs(sum) < 1e-10");
}}else {double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);

#line 194
{
  coord g;
  if (_attribute[s.i].gradient)
    {
      g.x = _attribute[s.i].gradient (val(s,-1,0,0), val(s,0,0,0), val(s,1,0,0));
      
#line 198
g.y = _attribute[s.i].gradient (val(s,0,-1,0), val(s,0,0,0), val(s,0,1,0));}
  else
    {
      g.x = (val(s,1,0,0) - val(s,-1,0,0))/2.;
      
#line 201
g.y = (val(s,0,1,0) - val(s,0,-1,0))/2.;}

  double sc = val(s,0,0,0), cmc = 4.*_const_cm, sum = _const_cm*(1 << 2);
  {foreach_child() {
    val(s,0,0,0) = sc;
    
      val(s,0,0,0) += child.x*g.x*_const_cm/cmc;
      
#line 207
val(s,0,0,0) += child.y*g.y*_const_cm/cmc;
    sum -= _const_cm;
  }end_foreach_child()}
  if (!(fabs(sum) < 1e-10)) qassert ("/u/basilisk/src/grid/multigrid-common.h", 0, "fabs(sum) < 1e-10");
}}

#line 211
}

static inline void refine_reset (Point point, scalar v)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  {foreach_child()
    val(v,0,0,0) = 0.;end_foreach_child()}
}

static inline void refine_injection (Point point, scalar v)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  double val = val(v,0,0,0);
  {foreach_child()
    val(v,0,0,0) = val;end_foreach_child()}
}

static scalar multigrid_init_scalar (scalar s, const char * name)
{
  s = cartesian_init_scalar (s, name);
  _attribute[s.i].prolongation = refine_bilinear;
  _attribute[s.i].restriction = restriction_average;
  return s;
}

static scalar multigrid_init_vertex_scalar (scalar s, const char * name)
{
  s = cartesian_init_vertex_scalar (s, name);
  _attribute[s.i].restriction = restriction_vertex;
  return s;
}

static vector multigrid_init_face_vector (vector v, const char * name)
{
  v = cartesian_init_face_vector (v, name);
  
    _attribute[v.y.i].restriction = no_restriction;
    
#line 245
_attribute[v.x.i].restriction = no_restriction;
  _attribute[v.x.i].restriction = restriction_face;
  return v;
}

void multigrid_debug (Point point)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  cartesian_debug (point);

  FILE * plot = fopen ("plot", "a");
  if (point.level > 0) {
    char name[80] = "coarse";
    if (pid() > 0)
      sprintf (name, "coarse-%d", pid());
    FILE * fp = fopen (name, "w");
#line 271 "/u/basilisk/src/grid/multigrid-common.h"
      double xc = x - child.x*Delta/2., yc = y - child.y*Delta/2.;
      for (int k = 0; k <= 1; k++)
 for (int l = 0; l <= 1; l++) {
   {scalar*_i=(scalar*)( all);if(_i)for(scalar v=*_i;(&v)->i>=0;v=*++_i){
     fprintf (fp, "%g %g %g ",
       xc + k*child.x*Delta*2. + _attribute[v.i].d.x*Delta,
       yc + l*child.y*Delta*2. + _attribute[v.i].d.y*Delta,
       coarse(v,k*child.x,l*child.y,0));}}
   fputc ('\n', fp);
 }
      fprintf (ferr, ", '%s' u 1+3*v:2+3*v:3+3*v w labels tc lt 3 t ''", name);
      fprintf (plot, ", '%s' u 1+3*v:2+3*v:3+3*v w labels tc lt 3 t ''", name);
#line 302 "/u/basilisk/src/grid/multigrid-common.h"
    fclose (fp);
  }

  if (is_coarse()) {
    char name[80] = "fine";
    if (pid() > 0)
      sprintf (name, "fine-%d", pid());
    FILE * fp = fopen (name, "w");
#line 324 "/u/basilisk/src/grid/multigrid-common.h"
      double xf = x - Delta/4., yf = y - Delta/4.;
      for (int k = -2; k <= 3; k++)
 for (int l = -2; l <= 3; l++) {
   {scalar*_i=(scalar*)( all);if(_i)for(scalar v=*_i;(&v)->i>=0;v=*++_i){ {
     fprintf (fp, "%g %g ",
       xf + k*Delta/2. + _attribute[v.i].d.x*Delta/4.,
       yf + l*Delta/2. + _attribute[v.i].d.y*Delta/4.);
     if (allocated_child(k,l,0))
       fprintf (fp, "%g ", fine(v,k,l,0));
     else
       fputs ("n/a ", fp);
   }}}
   fputc ('\n', fp);
 }
      fprintf (ferr, ", '%s' u 1+3*v:2+3*v:3+3*v w labels tc lt 2 t ''", name);
      fprintf (plot, ", '%s' u 1+3*v:2+3*v:3+3*v w labels tc lt 2 t ''", name);
#line 362 "/u/basilisk/src/grid/multigrid-common.h"
    fclose (fp);
  }
  fflush (ferr);
  fclose (plot);
}

static void multigrid_restriction (scalar * list)
{
  scalar * listdef = NULL, * listc = NULL, * list2 = NULL;
  {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    if (!is_constant (s) && _attribute[s.i].block > 0) {
      if (_attribute[s.i].restriction == restriction_average) {
 listdef = list_add (listdef, s);
 list2 = list_add (list2, s);
      }
      else if (_attribute[s.i].restriction != no_restriction) {
 listc = list_add (listc, s);
 if (_attribute[s.i].face)
   {
     list2 = list_add (list2, _attribute[s.i].v.x);
     
#line 381
list2 = list_add (list2, _attribute[s.i].v.y);}
 else
   list2 = list_add (list2, s);
      }
    }}}

  if (listdef || listc) {
    for (int l = depth() - 1; l >= 0; l--) {
      {foreach_coarse_level(l) {
 {scalar*_i=(scalar*)( listdef);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
  
     restriction_average (point, s);}}
 {scalar*_i=(scalar*)( listc);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {
  
     _attribute[s.i].restriction (point, s);
 }}}
      }end_foreach_coarse_level();}
      { Boundary ** _i = boundaries, * _b; while (_i && (_b = *_i++)) if (_b->level) _b->level (_b, list2, l); };
    }
    pfree (listdef,__func__,__FILE__,0);
    pfree (listc,__func__,__FILE__,0);
    pfree (list2,__func__,__FILE__,0);
  }
}

void multigrid_methods()
{
  cartesian_methods();
  debug = multigrid_debug;
  init_scalar = multigrid_init_scalar;
  init_vertex_scalar = multigrid_init_vertex_scalar;
  init_face_vector = multigrid_init_face_vector;
  restriction = multigrid_restriction;
}







void subtree_size (scalar size, bool leaves)
{




  foreach_stencil()
    {_stencil_val_a(size,0,0,0);  }end_foreach_stencil();




  {
#line 428
foreach()
    val(size,0,0,0) = 1;end_foreach();}





  { Boundary ** _i = boundaries, * _b; while (_i && (_b = *_i++)) if (_b->restriction) _b->restriction (_b, ((scalar[]){size,{-1}}), depth()); };
  for (int l = depth() - 1; l >= 0; l--) {
    {foreach_coarse_level(l) {
      double sum = !leaves;
      {foreach_child()
 sum += val(size,0,0,0);end_foreach_child()}
      val(size,0,0,0) = sum;
    }end_foreach_coarse_level();}
    { Boundary ** _i = boundaries, * _b; while (_i && (_b = *_i++)) if (_b->restriction) _b->restriction (_b, ((scalar[]){size,{-1}}), l); };
  }
}
#line 5 "/u/basilisk/src/grid/tree-common.h"




#line 21 "/u/basilisk/src/grid/tree-common.h"
int refine_cell (Point point, scalar * list, int flag, Cache * refined)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  int nr = 0;


  if (level > 0)
    for (int k = 0; k != 2*child.x; k += child.x)

      for (int l = 0; l != 2*child.y; l += child.y)




   if (aparent(k,l,m).pid >= 0 && is_leaf(aparent(k,l,m))) {
     Point p = point;


     p.level = point.level - 1;
     p.i = (point.i + 2)/2 + k;
     do { if (p.i < 2) p.i += 1 << p.level; else if (p.i >= 2 + (1 << p.level)) p.i -= 1 << p.level; } while(0);

       p.j = (point.j + 2)/2 + l;
       do { if (p.j < 2) p.j += 1 << p.level; else if (p.j >= 2 + (1 << p.level)) p.j -= 1 << p.level; } while(0);





     nr += refine_cell (p, list, flag, refined);
     aparent(k,l,m).flags |= flag;
   }



  increment_neighbors (point);

  int cflag = is_active(cell) ? (active|leaf) : leaf;
  {foreach_child()
    cell.flags |= cflag;end_foreach_child()}


  {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    if (is_local(cell) || _attribute[s.i].face)
      _attribute[s.i].refine (point, s);}}


  cell.flags &= ~leaf;

#if _MPI
  if (is_border(cell)) {
    {foreach_child() {
      bool bord = false;
      {foreach_neighbor() {
 if (!is_local(cell) || (level > 0 && !is_local(aparent(0,0,0)))) {
   bord = true; foreach_neighbor_break();
 }
 if (is_refined_check())
   {foreach_child()
     if (!is_local(cell)) {
       bord = true; foreach_child_break();
     }end_foreach_child()}
 if (bord)
   foreach_neighbor_break();
      }end_foreach_neighbor()}
      if (bord)
 cell.flags |= border;
    }end_foreach_child()}
    if (refined)
      cache_append (refined, point, cell.flags);
    nr++;
  }
#endif
  return nr;
}





bool coarsen_cell (Point point, scalar * list)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;



  int pid = cell.pid;
  {foreach_child()
    if (cell.neighbors || (cell.pid < 0 && cell.pid != pid))
      return false;end_foreach_child()}



  {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {
    _attribute[s.i].restriction (point, s);
    if (_attribute[s.i].coarsen)
      _attribute[s.i].coarsen (point, s);
  }}}


  cell.flags |= leaf;


  decrement_neighbors (point);

#if _MPI
  if (!is_local(cell)) {
    cell.flags &= ~(active|border);
    {foreach_neighbor(1)
      if (cell.neighbors)
 {foreach_child()
   if (allocated(0,0,0) && is_local(cell) && !is_border(cell))
     cell.flags |= border;end_foreach_child()}end_foreach_neighbor()}
  }
#endif

  return true;
}

void coarsen_cell_recursive (Point point, scalar * list)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;


  {foreach_child()
    if (cell.neighbors)
      {foreach_neighbor(1)
 if ((!is_leaf (cell) && cell.neighbors && cell.pid >= 0))
   coarsen_cell_recursive (point, list);end_foreach_neighbor()}end_foreach_child()}

  if (!(coarsen_cell (point, list))) qassert ("/u/basilisk/src/grid/tree-common.h", 0, "coarsen_cell (point, list)");
}

void mpi_boundary_refine (scalar *);
void mpi_boundary_coarsen (int, int);
void mpi_boundary_update (scalar *);

typedef struct {
  int nc, nf;
} astats;

struct Adapt {
  scalar * slist;
  double * max;
  int maxlevel;
  int minlevel;
  scalar * list;
};

     
astats adapt_wavelet (struct Adapt p)
{tracing("adapt_wavelet","/u/basilisk/src/grid/tree-common.h",0);
  scalar * list = p.list;

  if (is_constant(cm)) {
    if (list == NULL || list == all)
      list = list_copy (all);
    boundary_internal ((scalar *)list, "/u/basilisk/src/grid/tree-common.h", 0);
    restriction (p.slist);
  }
  else {
    if (list == NULL || list == all) {
      list = list_copy (((scalar[]){cm, fm.x,fm.y,{-1}}));
      {scalar*_i=(scalar*)( all);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
 list = list_add (list, s);}}
    }
    boundary_internal ((scalar *)list, "/u/basilisk/src/grid/tree-common.h", 0);
    scalar * listr = list_concat (p.slist, ((scalar[]){cm,{-1}}));
    restriction (listr);
    pfree (listr,__func__,__FILE__,0);
  }

  astats st = {0, 0};
  scalar * listc = NULL;
  {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    if (!is_constant(s) && _attribute[s.i].restriction != no_restriction)
      listc = list_add (listc, s);}}


  if (p.minlevel < 1)
    p.minlevel = 1;
  ((Tree *)grid)->refined.n = 0;
  static const int refined = 1 << user, too_fine = 1 << (user + 1);
  {foreach_cell() {
    if (is_active(cell)) {
      static const int too_coarse = 1 << (user + 2);
      if (is_leaf (cell)) {
 if (cell.flags & too_coarse) {
   cell.flags &= ~too_coarse;
   refine_cell (point, listc, refined, &((Tree *)grid)->refined);
   st.nf++;
 }
 continue;
      }
      else {
 if (cell.flags & refined) {

   cell.flags &= ~too_coarse;
   continue;
 }

 bool local = is_local(cell);
 if (!local)
   {foreach_child()
     if (is_local(cell)) {
       local = true; foreach_child_break();
     }end_foreach_child()}
 if (local) {
   int i = 0;
   static const int just_fine = 1 << (user + 3);
   {scalar*_i=(scalar*)( p.slist);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {
     double max = p.max[i++], sc[1 << 2];
     int c = 0;
     {foreach_child()
       sc[c++] = val(s,0,0,0);end_foreach_child()}
     _attribute[s.i].prolongation (point, s);
     c = 0;
     {foreach_child() {
       double e = fabs(sc[c] - val(s,0,0,0));
       if (e > max && level < p.maxlevel) {
  cell.flags &= ~too_fine;
  cell.flags |= too_coarse;
       }
       else if ((e <= max/1.5 || level > p.maxlevel) &&
         !(cell.flags & (too_coarse|just_fine))) {
  if (level >= p.minlevel)
    cell.flags |= too_fine;
       }
       else if (!(cell.flags & too_coarse)) {
  cell.flags &= ~too_fine;
  cell.flags |= just_fine;
       }
       val(s,0,0,0) = sc[c++];
     }end_foreach_child()}
   }}}
   {foreach_child() {
     cell.flags &= ~just_fine;
     if (!is_leaf(cell)) {
       cell.flags &= ~too_coarse;
       if (level >= p.maxlevel)
  cell.flags |= too_fine;
     }
     else if (!is_active(cell))
       cell.flags &= ~too_coarse;
   }end_foreach_child()}
 }
      }
    }
    else
      continue;
  }end_foreach_cell();}
  mpi_boundary_refine (listc);



  for (int l = depth(); l >= 0; l--) {
    {foreach_cell()
      if (!(cell.pid < 0)) {
 if (level == l) {
   if (!is_leaf(cell)) {
     if (cell.flags & refined)

       cell.flags &= ~(refined|too_fine);
     else if (cell.flags & too_fine) {
       if (is_local(cell) && coarsen_cell (point, listc))
  st.nc++;
       cell.flags &= ~too_fine;
     }
   }
   if (cell.flags & too_fine)
     cell.flags &= ~too_fine;
   else if (level > 0 && (aparent(0,0,0).flags & too_fine))
     aparent(0,0,0).flags &= ~too_fine;
   continue;
 }
 else if (is_leaf(cell))
   continue;
      }end_foreach_cell();}
    mpi_boundary_coarsen (l, too_fine);
  }
  pfree (listc,__func__,__FILE__,0);

  mpi_all_reduce (st.nf, MPI_INT, MPI_SUM);
  mpi_all_reduce (st.nc, MPI_INT, MPI_SUM);
  if (st.nc || st.nf)
    mpi_boundary_update (list);

  if (list != p.list)
    pfree (list,__func__,__FILE__,0);

  {end_tracing("adapt_wavelet","/u/basilisk/src/grid/tree-common.h",0);return st;}
end_tracing("adapt_wavelet","/u/basilisk/src/grid/tree-common.h",0);}
#line 331 "/u/basilisk/src/grid/tree-common.h"
static void refine_level (int depth)
{
  int refined;
  do {
    refined = 0;
    ((Tree *)grid)->refined.n = 0;
    {foreach_leaf()
      if (level < depth) {
 refine_cell (point, NULL, 0, &((Tree *)grid)->refined);
 refined++;
 continue;
      }end_foreach_leaf();}
    mpi_all_reduce (refined, MPI_INT, MPI_SUM);
    if (refined) {
      mpi_boundary_refine (NULL);
      mpi_boundary_update (NULL);
    }
  } while (refined);
}
#line 376 "/u/basilisk/src/grid/tree-common.h"
     
static void halo_face (vectorl vl)
{tracing("halo_face","/u/basilisk/src/grid/tree-common.h",0);
  
    {scalar*_i=(scalar*)( vl.x);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      _attribute[s.i].dirty = 2;}}
    
#line 380
{scalar*_i=(scalar*)( vl.y);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      _attribute[s.i].dirty = 2;}}

  for (int l = depth() - 1; l >= 0; l--)
    {foreach_halo (prolongation, l)
      {
        if (vl.x) {
#line 395 "/u/basilisk/src/grid/tree-common.h"
   if ((!is_leaf (neighbor(-1,0,0)) && neighbor(-1,0,0).neighbors && neighbor(-1,0,0).pid >= 0))
     {scalar*_i=(scalar*)( vl.x);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
       val(s,0,0,0) = (fine(s,0,0,0) + fine(s,0,1,0))/2.;}}
   if ((!is_leaf (neighbor(1,0,0)) && neighbor(1,0,0).neighbors && neighbor(1,0,0).pid >= 0))
     {scalar*_i=(scalar*)( vl.x);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
       val(s,1,0,0) = (fine(s,2,0,0) + fine(s,2,1,0))/2.;}}
#line 411 "/u/basilisk/src/grid/tree-common.h"
 }
        
#line 386
if (vl.y) {
#line 395 "/u/basilisk/src/grid/tree-common.h"
   if ((!is_leaf (neighbor(0,-1,0)) && neighbor(0,-1,0).neighbors && neighbor(0,-1,0).pid >= 0))
     {scalar*_i=(scalar*)( vl.y);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
       val(s,0,0,0) = (fine(s,0,0,0) + fine(s,1,0,0))/2.;}}
   if ((!is_leaf (neighbor(0,1,0)) && neighbor(0,1,0).neighbors && neighbor(0,1,0).pid >= 0))
     {scalar*_i=(scalar*)( vl.y);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
       val(s,0,1,0) = (fine(s,0,2,0) + fine(s,1,2,0))/2.;}}
#line 411 "/u/basilisk/src/grid/tree-common.h"
 }}end_foreach_halo();}
end_tracing("halo_face","/u/basilisk/src/grid/tree-common.h",0);}



static scalar tree_init_scalar (scalar s, const char * name)
{
  s = multigrid_init_scalar (s, name);
  _attribute[s.i].refine = _attribute[s.i].prolongation;
  return s;
}

static void prolongation_vertex (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;

  fine(s,1,1,0) = (val(s,0,0,0) + val(s,1,0,0) + val(s,0,1,0) + val(s,1,1,0))/4.;





  for (int i = 0; i <= 1; i++) {
    for (int j = 0; j <= 1; j++)





      if (allocated_child(2*i,2*j,0))
 fine(s,2*i,2*j,0) = val(s,i,j,0);


    
      if (neighbor(i,0,0).neighbors) {

 fine(s,2*i,1,0) = (val(s,i,0,0) + val(s,i,1,0))/2.;
#line 456 "/u/basilisk/src/grid/tree-common.h"
      }
      
#line 444
if (neighbor(0,i,0).neighbors) {

 fine(s,1,2*i,0) = (val(s,0,i,0) + val(s,1,i,0))/2.;
#line 456 "/u/basilisk/src/grid/tree-common.h"
      }
  }
}

static scalar tree_init_vertex_scalar (scalar s, const char * name)
{
  s = multigrid_init_vertex_scalar (s, name);
  _attribute[s.i].refine = _attribute[s.i].prolongation = prolongation_vertex;
  return s;
}


static void refine_face_x (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  vector v = _attribute[s.i].v;

  if (!(!is_leaf (neighbor(-1,0,0)) && neighbor(-1,0,0).neighbors && neighbor(-1,0,0).pid >= 0) &&
      (is_local(cell) || is_local(neighbor(-1,0,0)))) {
    double g1 = (val(v.x,0,+1,0) - val(v.x,0,-1,0))/8.;
    for (int j = 0; j <= 1; j++)
      fine(v.x,0,j,0) = val(v.x,0,0,0) + (2*j - 1)*g1;
  }
  if (!(!is_leaf (neighbor(1,0,0)) && neighbor(1,0,0).neighbors && neighbor(1,0,0).pid >= 0) && neighbor(1,0,0).neighbors &&
      (is_local(cell) || is_local(neighbor(1,0,0)))) {
    double g1 = (val(v.x,1,+1,0) - val(v.x,1,-1,0))/8.;
    for (int j = 0; j <= 1; j++)
      fine(v.x,2,j,0) = val(v.x,1,0,0) + (2*j - 1)*g1;
  }
  if (is_local(cell)) {
    double g1 = (val(v.x,0,+1,0) - val(v.x,0,-1,0) + val(v.x,1,+1,0) - val(v.x,1,-1,0))/16.;
    for (int j = 0; j <= 1; j++)
      fine(v.x,1,j,0) = (val(v.x,0,0,0) + val(v.x,1,0,0))/2. + (2*j - 1)*g1;
  }
#line 514 "/u/basilisk/src/grid/tree-common.h"
}

#line 468
static void refine_face_y (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  vector v = _attribute[s.i].v;

  if (!(!is_leaf (neighbor(0,-1,0)) && neighbor(0,-1,0).neighbors && neighbor(0,-1,0).pid >= 0) &&
      (is_local(cell) || is_local(neighbor(0,-1,0)))) {
    double g1 = (val(v.y,+1,0,0) - val(v.y,-1,0,0))/8.;
    for (int j = 0; j <= 1; j++)
      fine(v.y,j,0,0) = val(v.y,0,0,0) + (2*j - 1)*g1;
  }
  if (!(!is_leaf (neighbor(0,1,0)) && neighbor(0,1,0).neighbors && neighbor(0,1,0).pid >= 0) && neighbor(0,1,0).neighbors &&
      (is_local(cell) || is_local(neighbor(0,1,0)))) {
    double g1 = (val(v.y,+1,1,0) - val(v.y,-1,1,0))/8.;
    for (int j = 0; j <= 1; j++)
      fine(v.y,j,2,0) = val(v.y,0,1,0) + (2*j - 1)*g1;
  }
  if (is_local(cell)) {
    double g1 = (val(v.y,+1,0,0) - val(v.y,-1,0,0) + val(v.y,+1,1,0) - val(v.y,-1,1,0))/16.;
    for (int j = 0; j <= 1; j++)
      fine(v.y,j,1,0) = (val(v.y,0,0,0) + val(v.y,0,1,0))/2. + (2*j - 1)*g1;
  }
#line 514 "/u/basilisk/src/grid/tree-common.h"
}

void refine_face (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  vector v = _attribute[s.i].v;
  
    _attribute[v.x.i].prolongation (point, v.x);
    
#line 520
_attribute[v.y.i].prolongation (point, v.y);
}

void refine_face_solenoidal (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  refine_face (point, s);

  if (is_local(cell)) {

    vector v = _attribute[s.i].v;
    double d[1 << 2], p[1 << 2];
    int i = 0;
    {foreach_child() {
      d[i] = 0.;
      
 d[i] += val(v.x,1,0,0) - val(v.x,0,0,0);
 
#line 535
d[i] += val(v.y,0,1,0) - val(v.y,0,0,0);
      i++;
    }end_foreach_child()}

    p[0] = 0.;
    p[1] = (3.*d[3] + d[0])/4. + d[2]/2.;
    p[2] = (d[3] + 3.*d[0])/4. + d[2]/2.;
    p[3] = (d[3] + d[0])/2. + d[2];
    fine(v.x,1,1,0) += p[1] - p[0];
    fine(v.x,1,0,0) += p[3] - p[2];
    fine(v.y,0,1,0) += p[0] - p[2];
    fine(v.y,1,1,0) += p[1] - p[3];
#line 574 "/u/basilisk/src/grid/tree-common.h"
  }

}

vector tree_init_face_vector (vector v, const char * name)
{
  v = cartesian_init_face_vector (v, name);
  
    _attribute[v.x.i].restriction = _attribute[v.x.i].refine = no_restriction;
    
#line 582
_attribute[v.y.i].restriction = _attribute[v.y.i].refine = no_restriction;
  _attribute[v.x.i].restriction = restriction_face;
  _attribute[v.x.i].refine = refine_face;
  
    _attribute[v.x.i].prolongation = refine_face_x;
    
#line 586
_attribute[v.y.i].prolongation = refine_face_y;
  return v;
}

     
static void tree_boundary_level (scalar * list, int l)
{tracing("tree_boundary_level","/u/basilisk/src/grid/tree-common.h",0);
  int depth = l < 0 ? depth() : l;

  if (tree_is_full()) {
    { Boundary ** _i = boundaries, * _b; while (_i && (_b = *_i++)) if (_b->level) _b->level (_b, list, depth); };
    {end_tracing("tree_boundary_level","/u/basilisk/src/grid/tree-common.h",0);return;}
  }

  scalar * listdef = NULL, * listc = NULL, * list2 = NULL, * vlist = NULL;
  {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    if (!is_constant (s)) {
      if (_attribute[s.i].restriction == restriction_average) {
 listdef = list_add (listdef, s);
 list2 = list_add (list2, s);
      }
      else if (_attribute[s.i].restriction != no_restriction) {
 listc = list_add (listc, s);
 if (_attribute[s.i].face)
   {
     list2 = list_add (list2, _attribute[s.i].v.x);
     
#line 611
list2 = list_add (list2, _attribute[s.i].v.y);}
 else {
   list2 = list_add (list2, s);
   if (_attribute[s.i].restriction == restriction_vertex)
     vlist = list_add (vlist, s);
 }
      }
    }}}

  if (vlist)






    {foreach_vertex () {
      if ((!is_leaf (cell) && cell.neighbors && cell.pid >= 0) || (!is_leaf (neighbor(-1,0,0)) && neighbor(-1,0,0).neighbors && neighbor(-1,0,0).pid >= 0) ||
   (!is_leaf (neighbor(0,-1,0)) && neighbor(0,-1,0).neighbors && neighbor(0,-1,0).pid >= 0) || (!is_leaf (neighbor(-1,-1,0)) && neighbor(-1,-1,0).neighbors && neighbor(-1,-1,0).pid >= 0)) {

 {scalar*_i=(scalar*)( vlist);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
   val(s,0,0,0) = is_vertex (child(0,0,0)) ? fine(s,0,0,0) : HUGE;}}
      }
      else
 {
   if (child.y == 1 &&
       ((!is_leaf(cell) && !cell.neighbors && cell.pid >= 0) || (!is_leaf(neighbor(-1,0,0)) && !neighbor(-1,0,0).neighbors && neighbor(-1,0,0).pid >= 0))) {

     {scalar*_i=(scalar*)( vlist);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
       val(s,0,0,0) = is_vertex(neighbor(0,-1,0)) && is_vertex(neighbor(0,1,0)) ?
  (val(s,0,-1,0) + val(s,0,1,0))/2. : HUGE;}}
   }
   
#line 636
if (child.x == 1 &&
       ((!is_leaf(cell) && !cell.neighbors && cell.pid >= 0) || (!is_leaf(neighbor(0,-1,0)) && !neighbor(0,-1,0).neighbors && neighbor(0,-1,0).pid >= 0))) {

     {scalar*_i=(scalar*)( vlist);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
       val(s,0,0,0) = is_vertex(neighbor(-1,0,0)) && is_vertex(neighbor(1,0,0)) ?
  (val(s,-1,0,0) + val(s,1,0,0))/2. : HUGE;}}
   }}
    }end_foreach_vertex();}
#line 676 "/u/basilisk/src/grid/tree-common.h"
  pfree (vlist,__func__,__FILE__,0);

  if (listdef || listc) {
    { Boundary ** _i = boundaries, * _b; while (_i && (_b = *_i++)) if (_b->restriction) _b->restriction (_b, list2, depth); };
    for (int l = depth - 1; l >= 0; l--) {
      {foreach_coarse_level(l) {
 {scalar*_i=(scalar*)( listdef);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
   restriction_average (point, s);}}
 {scalar*_i=(scalar*)( listc);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
   _attribute[s.i].restriction (point, s);}}
      }end_foreach_coarse_level();}
      { Boundary ** _i = boundaries, * _b; while (_i && (_b = *_i++)) if (_b->restriction) _b->restriction (_b, list2, l); };
    }
    pfree (listdef,__func__,__FILE__,0);
    pfree (listc,__func__,__FILE__,0);
    pfree (list2,__func__,__FILE__,0);
  }

  scalar * listr = NULL;
  vector * listf = NULL;
  {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    if (!is_constant (s) && _attribute[s.i].refine != no_restriction) {
      if (_attribute[s.i].face)
 listf = vectors_add (listf, _attribute[s.i].v);
      else
 listr = list_add (listr, s);
    }}}

  if (listr || listf) {
    { Boundary ** _i = boundaries, * _b; while (_i && (_b = *_i++)) if (_b->level) _b->level (_b, list, 0); };
    for (int i = 0; i < depth; i++) {
      {foreach_halo (prolongation, i) {
 {scalar*_i=(scalar*)( listr);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
          _attribute[s.i].prolongation (point, s);}}
 {vector*_i=(vector*)( listf);if(_i)for(vector v=*_i;(&v)->x.i>=0;v=*++_i){
   {
     _attribute[v.x.i].prolongation (point, v.x);
     
#line 712
_attribute[v.y.i].prolongation (point, v.y);}}}
      }end_foreach_halo();}
      { Boundary ** _i = boundaries, * _b; while (_i && (_b = *_i++)) if (_b->level) _b->level (_b, list, i + 1); };
    }
    pfree (listr,__func__,__FILE__,0);
    pfree (listf,__func__,__FILE__,0);
  }
end_tracing("tree_boundary_level","/u/basilisk/src/grid/tree-common.h",0);}

double treex (Point point) {int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  if (level == 0)
    return 0;

  double i = 2*child.x - child.y;
  if (i <= 1 && i >= -1) i = -i;




  return treex(parent) + i/(1 << 2*(level - 1));
}

double treey (Point point) {int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  if (level == 0)
    return 0;
  return treey(parent) + 4./(1 << 2*(level - 1));
}

void output_tree (FILE * fp)
{
  {foreach_cell()
    if (cell.neighbors)
      {foreach_child()
 if (is_local(cell))
   fprintf (fp, "%g %g\n%g %g\n\n",
     treex(parent), treey(parent), treex(point), treey(point));end_foreach_child()}end_foreach_cell();}
}

     
void tree_check()
{tracing("tree_check","/u/basilisk/src/grid/tree-common.h",0);


  long nleaves = 0, nactive = 0;
  {foreach_cell_all() {
    if (is_leaf(cell)) {
      if (!(cell.pid >= 0)) qassert ("/u/basilisk/src/grid/tree-common.h", 0, "cell.pid >= 0");
      nleaves++;
    }
    if (is_local(cell))
      if (!(is_active(cell) || (!is_leaf(cell) && !cell.neighbors && cell.pid >= 0))) qassert ("/u/basilisk/src/grid/tree-common.h", 0, "is_active(cell) || is_prolongation(cell)");
    if (is_active(cell))
      nactive++;

    int neighbors = 0;
    {foreach_neighbor(1)
      if (allocated(0,0,0) && (!is_leaf (cell) && cell.neighbors && cell.pid >= 0))
 neighbors++;end_foreach_neighbor()}
    if (!(cell.neighbors == neighbors)) qassert ("/u/basilisk/src/grid/tree-common.h", 0, "cell.neighbors == neighbors");


    if (!cell.neighbors)
      if (!(!allocated_child(0,0,0))) qassert ("/u/basilisk/src/grid/tree-common.h", 0, "!allocated_child(0)");
  }end_foreach_cell_all();}


  long reachable = 0;
  {foreach_cell() {
    if (is_active(cell))
      reachable++;
    else
      continue;
  }end_foreach_cell();}
  if (!(nactive == reachable)) qassert ("/u/basilisk/src/grid/tree-common.h", 0, "nactive == reachable");


  reachable = 0;
  {foreach_cell()
    if (is_leaf(cell)) {
      reachable++;
      continue;
    }end_foreach_cell();}
  if (!(nleaves == reachable)) qassert ("/u/basilisk/src/grid/tree-common.h", 0, "nleaves == reachable");
end_tracing("tree_check","/u/basilisk/src/grid/tree-common.h",0);}

     
static void tree_restriction (scalar * list) {tracing("tree_restriction","/u/basilisk/src/grid/tree-common.h",0);
  boundary_internal ((scalar *)list, "/u/basilisk/src/grid/tree-common.h", 0);
  if (tree_is_full())
    multigrid_restriction (list);
end_tracing("tree_restriction","/u/basilisk/src/grid/tree-common.h",0);}

void tree_methods()
{
  multigrid_methods();
  init_scalar = tree_init_scalar;
  init_vertex_scalar = tree_init_vertex_scalar;
  init_face_vector = tree_init_face_vector;
  boundary_level = tree_boundary_level;
  boundary_face = halo_face;
  restriction = tree_restriction;
}
#line 1672 "/u/basilisk/src/grid/tree.h"


void tree_periodic (int dir)
{
  int depth = grid ? depth() : -1;
  if (grid)
    free_grid();
  periodic (dir);
  if (depth >= 0)
    init_grid (1 << depth);
}


#if _MPI
#line 1 "grid/tree-mpi.h"
#line 1 "/u/basilisk/src/grid/tree-mpi.h"

int debug_iteration = -1;

void debug_mpi (FILE * fp1);

typedef struct {
  CacheLevel * halo;
  void * buf;
  MPI_Request r;
  int depth;
  int pid;
  int maxdepth;
} Rcv;

typedef struct {
  Rcv * rcv;
  char * name;
  int npid;
} RcvPid;

typedef struct {
  RcvPid * rcv, * snd;
} SndRcv;

typedef struct {
  Boundary parent;

  SndRcv mpi_level, mpi_level_root, restriction;
  Array * send, * receive;
} MpiBoundary;

static void cache_level_init (CacheLevel * c)
{
  c->p = NULL;
  c->n = c->nm = 0;
}

static void rcv_append (Point point, Rcv * rcv)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  if (level > rcv->depth) {
    rcv->halo = (CacheLevel *) prealloc (rcv->halo, (level + 1)*sizeof(CacheLevel),__func__,__FILE__,0);
    for (int j = rcv->depth + 1; j <= level; j++)
      cache_level_init (&rcv->halo[j]);
    rcv->depth = level;
  }
  cache_level_append (&rcv->halo[level], point);
  if (level > rcv->maxdepth)
    rcv->maxdepth = level;
}

void rcv_print (Rcv * rcv, FILE * fp, const char * prefix)
{
  for (int l = 0; l <= rcv->depth; l++)
    if (rcv->halo[l].n > 0)
      {foreach_cache_level(rcv->halo[l], l)
 fprintf (fp, "%s%g %g %g %d %d\n", prefix, x, y, z, rcv->pid, level);end_foreach_cache_level();}
}

static void rcv_free_buf (Rcv * rcv)
{
  if (rcv->buf) {
    prof_start ("rcv_pid_receive");
    MPI_Wait (&rcv->r, MPI_STATUS_IGNORE);
    pfree (rcv->buf,__func__,__FILE__,0);
    rcv->buf = NULL;
    prof_stop();
  }
}

static void rcv_destroy (Rcv * rcv)
{
  rcv_free_buf (rcv);
  for (int i = 0; i <= rcv->depth; i++)
    if (rcv->halo[i].n > 0)
      pfree (rcv->halo[i].p,__func__,__FILE__,0);
  pfree (rcv->halo,__func__,__FILE__,0);
}

static RcvPid * rcv_pid_new (const char * name)
{
  RcvPid * r = ((RcvPid *) pcalloc (1, sizeof(RcvPid),__func__,__FILE__,0));
  r->name = pstrdup (name,__func__,__FILE__,0);
  return r;
}

static Rcv * rcv_pid_pointer (RcvPid * p, int pid)
{
  if (!(pid >= 0 && pid < npe())) qassert ("/u/basilisk/src/grid/tree-mpi.h", 0, "pid >= 0 && pid < npe()");

  int i;
  for (i = 0; i < p->npid; i++)
    if (pid == p->rcv[i].pid)
      break;

  if (i == p->npid) {
    p->rcv = (Rcv *) prealloc (p->rcv, (++p->npid)*sizeof(Rcv),__func__,__FILE__,0);
    Rcv * rcv = &p->rcv[p->npid-1];
    rcv->pid = pid;
    rcv->depth = rcv->maxdepth = 0;
    rcv->halo = ((CacheLevel *) pmalloc ((1)*sizeof(CacheLevel),__func__,__FILE__,0));
    rcv->buf = NULL;
    cache_level_init (&rcv->halo[0]);
  }
  return &p->rcv[i];
}

static void rcv_pid_append (RcvPid * p, int pid, Point point)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  rcv_append (point, rcv_pid_pointer (p, pid));
}

static void rcv_pid_append_pids (RcvPid * p, Array * pids)
{

  for (int i = 0; i < p->npid; i++) {
    int pid = p->rcv[i].pid, j, * a;
    for (j = 0, a = pids->p; j < pids->len/sizeof(int); j++,a++)
      if (*a == pid)
 break;
    if (j == pids->len/sizeof(int))
      array_append (pids, &pid, sizeof(int));
  }
}

void rcv_pid_write (RcvPid * p, const char * name)
{
  for (int i = 0; i < p->npid; i++) {
    Rcv * rcv = &p->rcv[i];
    char fname[80];
    sprintf (fname, "%s-%d-%d", name, pid(), rcv->pid);
    FILE * fp = fopen (fname, "w");
    rcv_print (rcv, fp, "");
    fclose (fp);
  }
}

static void rcv_pid_print (RcvPid * p, FILE * fp, const char * prefix)
{
  for (int i = 0; i < p->npid; i++)
    rcv_print (&p->rcv[i], fp, prefix);
}

static void rcv_pid_destroy (RcvPid * p)
{
  for (int i = 0; i < p->npid; i++)
    rcv_destroy (&p->rcv[i]);
  pfree (p->rcv,__func__,__FILE__,0);
  pfree (p->name,__func__,__FILE__,0);
  pfree (p,__func__,__FILE__,0);
}

static Boundary * mpi_boundary = NULL;






void debug_mpi (FILE * fp1);

static void apply_bc (Rcv * rcv, scalar * list, scalar * listv,
        vector * listf, int l, MPI_Status s)
{
  double * b = rcv->buf;
  {foreach_cache_level(rcv->halo[l], l) {
    {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {
      memcpy (&val(s,0,0,0), b, sizeof(double)*_attribute[s.i].block);
      b += _attribute[s.i].block;
    }}}
    {vector*_i=(vector*)( listf);if(_i)for(vector v=*_i;(&v)->x.i>=0;v=*++_i){
      { {
 memcpy (&val(v.x,0,0,0), b, sizeof(double)*_attribute[v.x.i].block);
 b += _attribute[v.x.i].block;
 if (*b != HUGE && allocated(1,0,0))
   memcpy (&val(v.x,1,0,0), b, sizeof(double)*_attribute[v.x.i].block);
 b += _attribute[v.x.i].block;
      } 
#line 171
{
 memcpy (&val(v.y,0,0,0), b, sizeof(double)*_attribute[v.y.i].block);
 b += _attribute[v.y.i].block;
 if (*b != HUGE && allocated(0,1,0))
   memcpy (&val(v.y,0,1,0), b, sizeof(double)*_attribute[v.y.i].block);
 b += _attribute[v.y.i].block;
      }}}}
    {scalar*_i=(scalar*)( listv);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {
      for (int i = 0; i <= 1; i++)
 for (int j = 0; j <= 1; j++)







          {
     if (*b != HUGE && allocated(i,j,0))
       memcpy (&val(s,i,j,0), b, sizeof(double)*_attribute[s.i].block);
     b += _attribute[s.i].block;
          }

    }}}
  }end_foreach_cache_level();}
  size_t size = b - (double *) rcv->buf;
  pfree (rcv->buf,__func__,__FILE__,0);
  rcv->buf = NULL;

  int rlen;
  MPI_Get_count (&s, MPI_DOUBLE, &rlen);
  if (rlen != size) {
    fprintf (ferr,
      "rlen (%d) != size (%ld), %d receiving from %d at level %d\n"
      "Calling debug_mpi(NULL)...\n"
      "Aborting...\n",
      rlen, size, pid(), rcv->pid, l);
    fflush (ferr);
    debug_mpi (NULL);
    MPI_Abort (MPI_COMM_WORLD, -2);
  }
}
#line 234 "/u/basilisk/src/grid/tree-mpi.h"
static void mpi_recv_check (void * buf, int count, MPI_Datatype datatype,
       int source, int tag,
       MPI_Comm comm, MPI_Status * status,
       const char * name)
{
#line 269 "/u/basilisk/src/grid/tree-mpi.h"
  int errorcode = MPI_Recv (buf, count, datatype, source, tag, comm, status);
  if (errorcode != MPI_SUCCESS) {
    char string[MPI_MAX_ERROR_STRING];
    int resultlen;
    MPI_Error_string (errorcode, string, &resultlen);
    fprintf (ferr,
      "ERROR MPI_Recv \"%s\" (count = %d, source = %d, tag = %d):\n%s\n"
      "Calling debug_mpi(NULL)...\n"
      "Aborting...\n",
      name, count, source, tag, string);
    fflush (ferr);
    debug_mpi (NULL);
    MPI_Abort (MPI_COMM_WORLD, -1);
  }





}

     
static int mpi_waitany (int count, MPI_Request array_of_requests[], int *indx,
   MPI_Status *status)
{tracing("mpi_waitany","/u/basilisk/src/grid/tree-mpi.h",0);
  { int _ret= MPI_Waitany (count, array_of_requests, indx, status);end_tracing("mpi_waitany","/u/basilisk/src/grid/tree-mpi.h",0);return _ret;}
end_tracing("mpi_waitany","/u/basilisk/src/grid/tree-mpi.h",0);}

static int list_lenb (scalar * list) {
  int len = 0;
  {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    len += _attribute[s.i].block;}}
  return len;
}

static int vectors_lenb (vector * list) {
  int len = 0;
  {vector*_i=(vector*)( list);if(_i)for(vector v=*_i;(&v)->x.i>=0;v=*++_i){
    len += _attribute[v.x.i].block;}}
  return len;
}

static void rcv_pid_receive (RcvPid * m, scalar * list, scalar * listv,
        vector * listf, int l)
{
  if (m->npid == 0)
    return;

  prof_start ("rcv_pid_receive");

  int len = list_lenb (list) + 2*2*vectors_lenb (listf) +
    (1 << 2)*list_lenb (listv);

  MPI_Request r[m->npid];
  Rcv * rrcv[m->npid];
  int nr = 0;
  for (int i = 0; i < m->npid; i++) {
    Rcv * rcv = &m->rcv[i];
    if (l <= rcv->depth && rcv->halo[l].n > 0) {
      if (!(!rcv->buf)) qassert ("/u/basilisk/src/grid/tree-mpi.h", 0, "!rcv->buf");
      rcv->buf = pmalloc (sizeof (double)*rcv->halo[l].n*len,__func__,__FILE__,0);






      MPI_Irecv (rcv->buf, rcv->halo[l].n*len, MPI_DOUBLE, rcv->pid,
   (l), MPI_COMM_WORLD, &r[nr]);
      rrcv[nr++] = rcv;






    }
  }


  if (nr > 0) {
    int i;
    MPI_Status s;
    mpi_waitany (nr, r, &i, &s);
    while (i != MPI_UNDEFINED) {
      Rcv * rcv = rrcv[i];
      if (!(l <= rcv->depth && rcv->halo[l].n > 0)) qassert ("/u/basilisk/src/grid/tree-mpi.h", 0, "l <= rcv->depth && rcv->halo[l].n > 0");
      if (!(rcv->buf)) qassert ("/u/basilisk/src/grid/tree-mpi.h", 0, "rcv->buf");
      apply_bc (rcv, list, listv, listf, l, s);
      mpi_waitany (nr, r, &i, &s);
    }
  }

  prof_stop();
}

     
static void rcv_pid_wait (RcvPid * m)
{tracing("rcv_pid_wait","/u/basilisk/src/grid/tree-mpi.h",0);

  for (int i = 0; i < m->npid; i++)
    rcv_free_buf (&m->rcv[i]);
end_tracing("rcv_pid_wait","/u/basilisk/src/grid/tree-mpi.h",0);}

static void rcv_pid_send (RcvPid * m, scalar * list, scalar * listv,
     vector * listf, int l)
{
  if (m->npid == 0)
    return;

  prof_start ("rcv_pid_send");

  int len = list_lenb (list) + 2*2*vectors_lenb (listf) +
    (1 << 2)*list_lenb (listv);


  for (int i = 0; i < m->npid; i++) {
    Rcv * rcv = &m->rcv[i];
    if (l <= rcv->depth && rcv->halo[l].n > 0) {
      if (!(!rcv->buf)) qassert ("/u/basilisk/src/grid/tree-mpi.h", 0, "!rcv->buf");
      rcv->buf = pmalloc (sizeof (double)*rcv->halo[l].n*len,__func__,__FILE__,0);
      double * b = rcv->buf;
      {foreach_cache_level(rcv->halo[l], l) {
 {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {
   memcpy (b, &val(s,0,0,0), sizeof(double)*_attribute[s.i].block);
   b += _attribute[s.i].block;
 }}}
 {vector*_i=(vector*)( listf);if(_i)for(vector v=*_i;(&v)->x.i>=0;v=*++_i){
   { {
     memcpy (b, &val(v.x,0,0,0), sizeof(double)*_attribute[v.x.i].block);
     b += _attribute[v.x.i].block;
     if (allocated(1,0,0))
       memcpy (b, &val(v.x,1,0,0), sizeof(double)*_attribute[v.x.i].block);
     else
       *b = HUGE;
     b += _attribute[v.x.i].block;
   } 
#line 397
{
     memcpy (b, &val(v.y,0,0,0), sizeof(double)*_attribute[v.y.i].block);
     b += _attribute[v.y.i].block;
     if (allocated(0,1,0))
       memcpy (b, &val(v.y,0,1,0), sizeof(double)*_attribute[v.y.i].block);
     else
       *b = HUGE;
     b += _attribute[v.y.i].block;
   }}}}
 {scalar*_i=(scalar*)( listv);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {
   for (int i = 0; i <= 1; i++)
     for (int j = 0; j <= 1; j++)
#line 418 "/u/basilisk/src/grid/tree-mpi.h"
       {
  if (allocated(i,j,0))
    memcpy (b, &val(s,i,j,0), sizeof(double)*_attribute[s.i].block);
  else
    *b = HUGE;
  b += _attribute[s.i].block;
       }

 }}}
      }end_foreach_cache_level();}





      MPI_Isend (rcv->buf, (b - (double *) rcv->buf),
   MPI_DOUBLE, rcv->pid, (l), MPI_COMM_WORLD,
   &rcv->r);
    }
  }

  prof_stop();
}

static void rcv_pid_sync (SndRcv * m, scalar * list, int l)
{
  scalar * listr = NULL, * listv = NULL;
  vector * listf = NULL;
  {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    if (!is_constant(s) && _attribute[s.i].block > 0) {
      if (_attribute[s.i].face)
 listf = vectors_add (listf, _attribute[s.i].v);
      else if (_attribute[s.i].restriction == restriction_vertex)
 listv = list_add (listv, s);
      else
 listr = list_add (listr, s);
    }}}
  rcv_pid_send (m->snd, listr, listv, listf, l);
  rcv_pid_receive (m->rcv, listr, listv, listf, l);
  rcv_pid_wait (m->snd);
  pfree (listr,__func__,__FILE__,0);
  pfree (listf,__func__,__FILE__,0);
  pfree (listv,__func__,__FILE__,0);
}

static void snd_rcv_destroy (SndRcv * m)
{
  rcv_pid_destroy (m->rcv);
  rcv_pid_destroy (m->snd);
}

static void snd_rcv_init (SndRcv * m, const char * name)
{
  char s[strlen(name) + 5];
  strcpy (s, name);
  strcat (s, ".rcv");
  m->rcv = rcv_pid_new (s);
  strcpy (s, name);
  strcat (s, ".snd");
  m->snd = rcv_pid_new (s);
}

static void mpi_boundary_destroy (Boundary * b)
{
  MpiBoundary * m = (MpiBoundary *) b;
  snd_rcv_destroy (&m->mpi_level);
  snd_rcv_destroy (&m->mpi_level_root);
  snd_rcv_destroy (&m->restriction);
  array_free (m->send);
  array_free (m->receive);
  pfree (m,__func__,__FILE__,0);
}

     
static void mpi_boundary_level (const Boundary * b, scalar * list, int l)
{tracing("mpi_boundary_level","/u/basilisk/src/grid/tree-mpi.h",0);
  MpiBoundary * m = (MpiBoundary *) b;
  rcv_pid_sync (&m->mpi_level, list, l);
  rcv_pid_sync (&m->mpi_level_root, list, l);
end_tracing("mpi_boundary_level","/u/basilisk/src/grid/tree-mpi.h",0);}

     
static void mpi_boundary_restriction (const Boundary * b, scalar * list, int l)
{tracing("mpi_boundary_restriction","/u/basilisk/src/grid/tree-mpi.h",0);
  MpiBoundary * m = (MpiBoundary *) b;
  rcv_pid_sync (&m->restriction, list, l);
end_tracing("mpi_boundary_restriction","/u/basilisk/src/grid/tree-mpi.h",0);}

void mpi_boundary_new()
{
  mpi_boundary = (Boundary *) ((MpiBoundary *) pcalloc (1, sizeof(MpiBoundary),__func__,__FILE__,0));
  mpi_boundary->destroy = mpi_boundary_destroy;
  mpi_boundary->level = mpi_boundary_level;
  mpi_boundary->restriction = mpi_boundary_restriction;
  MpiBoundary * mpi = (MpiBoundary *) mpi_boundary;
  snd_rcv_init (&mpi->mpi_level, "mpi_level");
  snd_rcv_init (&mpi->mpi_level_root, "mpi_level_root");
  snd_rcv_init (&mpi->restriction, "restriction");
  mpi->send = array_new();
  mpi->receive = array_new();
  add_boundary (mpi_boundary);
}

static FILE * fopen_prefix (FILE * fp, const char * name, char * prefix)
{
  if (fp) {
    sprintf (prefix, "%s-%d ", name, pid());
    return fp;
  }
  else {
    strcpy (prefix, "");
    char fname[80];
    if (debug_iteration >= 0)
      sprintf (fname, "%s-%d-%d", name, debug_iteration, pid());
    else
      sprintf (fname, "%s-%d", name, pid());
    return fopen (fname, "w");
  }
}

void debug_mpi (FILE * fp1)
{
  void output_cells_internal (FILE * fp);

  char prefix[80];
  FILE * fp;


  if (fp1 == NULL) {
    char name[80];
    sprintf (name, "halo-%d", pid()); remove (name);
    sprintf (name, "cells-%d", pid()); remove (name);
    sprintf (name, "faces-%d", pid()); remove (name);
    sprintf (name, "vertices-%d", pid()); remove (name);
    sprintf (name, "neighbors-%d", pid()); remove (name);
    sprintf (name, "mpi-level-rcv-%d", pid()); remove (name);
    sprintf (name, "mpi-level-snd-%d", pid()); remove (name);
    sprintf (name, "mpi-level-root-rcv-%d", pid()); remove (name);
    sprintf (name, "mpi-level-root-snd-%d", pid()); remove (name);
    sprintf (name, "mpi-restriction-rcv-%d", pid()); remove (name);
    sprintf (name, "mpi-restriction-snd-%d", pid()); remove (name);
    sprintf (name, "mpi-border-%d", pid()); remove (name);
    sprintf (name, "exterior-%d", pid()); remove (name);
    sprintf (name, "depth-%d", pid()); remove (name);
    sprintf (name, "refined-%d", pid()); remove (name);
  }


  fp = fopen_prefix (fp1, "halo", prefix);
  for (int l = 0; l < depth(); l++)
    {foreach_halo (prolongation, l)
      {foreach_child()
        fprintf (fp, "%s%g %g %g %d\n", prefix, x, y, z, level);end_foreach_child()}end_foreach_halo();}
  if (!fp1)
    fclose (fp);

  if (!fp1) {
    fp = fopen_prefix (fp1, "cells", prefix);
    output_cells_internal (fp);
    fclose (fp);
  }

  fp = fopen_prefix (fp1, "faces", prefix);
  {foreach_face_generic(){is_face_x(){
    fprintf (fp, "%s%g %g %g %d\n", prefix, x, y, z, level);}end_is_face_x()
#line 581
is_face_y(){
    fprintf (fp, "%s%g %g %g %d\n", prefix, x, y, z, level);}end_is_face_y()}end_foreach_face_generic();}
  if (!fp1)
    fclose (fp);

  fp = fopen_prefix (fp1, "vertices", prefix);
  {foreach_vertex()
    fprintf (fp, "%s%g %g %g %d\n", prefix, x, y, z, level);end_foreach_vertex();}
  if (!fp1)
    fclose (fp);

  fp = fopen_prefix (fp1, "neighbors", prefix);
  {foreach() {
    int n = 0;
    {foreach_neighbor(1)
      if ((!is_leaf (cell) && cell.neighbors && cell.pid >= 0))
 n++;end_foreach_neighbor()}
    fprintf (fp, "%s%g %g %g %d\n", prefix, x, y, z, cell.neighbors);
    if (!(cell.neighbors == n)) qassert ("/u/basilisk/src/grid/tree-mpi.h", 0, "cell.neighbors == n");
  }end_foreach();}
  if (!fp1)
    fclose (fp);

  MpiBoundary * mpi = (MpiBoundary *) mpi_boundary;

  fp = fopen_prefix (fp1, "mpi-level-rcv", prefix);
  rcv_pid_print (mpi->mpi_level.rcv, fp, prefix);
  if (!fp1)
    fclose (fp);

  fp = fopen_prefix (fp1, "mpi-level-root-rcv", prefix);
  rcv_pid_print (mpi->mpi_level_root.rcv, fp, prefix);
  if (!fp1)
    fclose (fp);

  fp = fopen_prefix (fp1, "mpi-restriction-rcv", prefix);
  rcv_pid_print (mpi->restriction.rcv, fp, prefix);
  if (!fp1)
    fclose (fp);

  fp = fopen_prefix (fp1, "mpi-level-snd", prefix);
  rcv_pid_print (mpi->mpi_level.snd, fp, prefix);
  if (!fp1)
    fclose (fp);

  fp = fopen_prefix (fp1, "mpi-level-root-snd", prefix);
  rcv_pid_print (mpi->mpi_level_root.snd, fp, prefix);
  if (!fp1)
    fclose (fp);

  fp = fopen_prefix (fp1, "mpi-restriction-snd", prefix);
  rcv_pid_print (mpi->restriction.snd, fp, prefix);
  if (!fp1)
    fclose (fp);

  fp = fopen_prefix (fp1, "mpi-border", prefix);
  {foreach_cell() {
    if (is_border(cell))
      fprintf (fp, "%s%g %g %g %d %d %d\n",
        prefix, x, y, z, level, cell.neighbors, cell.pid);
    else
      continue;
    if (is_leaf(cell))
      continue;
  }end_foreach_cell();}
  if (!fp1)
    fclose (fp);

  fp = fopen_prefix (fp1, "exterior", prefix);
  {foreach_cell() {
    if (!is_local(cell))
      fprintf (fp, "%s%g %g %g %d %d %d %d\n",
        prefix, x, y, z, level, cell.neighbors,
        cell.pid, cell.flags & leaf);






  }end_foreach_cell();}
  if (!fp1)
    fclose (fp);

  fp = fopen_prefix (fp1, "depth", prefix);
  fprintf (fp, "depth: %d %d\n", pid(), depth());
  fprintf (fp, "======= mpi_level.snd ======\n");
  RcvPid * snd = mpi->mpi_level.snd;
  for (int i = 0; i < snd->npid; i++)
    fprintf (fp, "%d %d %d\n", pid(), snd->rcv[i].pid, snd->rcv[i].maxdepth);
  fprintf (fp, "======= mpi_level.rcv ======\n");
  snd = mpi->mpi_level.rcv;
  for (int i = 0; i < snd->npid; i++)
    fprintf (fp, "%d %d %d\n", pid(), snd->rcv[i].pid, snd->rcv[i].maxdepth);
  if (!fp1)
    fclose (fp);

  fp = fopen_prefix (fp1, "refined", prefix);
  {foreach_cache (((Tree *)grid)->refined)
    fprintf (fp, "%s%g %g %g %d\n", prefix, x, y, z, level);end_foreach_cache();}
  if (!fp1)
    fclose (fp);
}

static void snd_rcv_free (SndRcv * p)
{
  char name[strlen(p->rcv->name) + 1];
  strcpy (name, p->rcv->name);
  rcv_pid_destroy (p->rcv);
  p->rcv = rcv_pid_new (name);
  strcpy (name, p->snd->name);
  rcv_pid_destroy (p->snd);
  p->snd = rcv_pid_new (name);
}

static bool is_root (Point point)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  if ((!is_leaf (cell) && cell.neighbors && cell.pid >= 0))
    {foreach_child()
      if (is_local(cell))
 return true;end_foreach_child()}
  return false;
}


static bool is_local_prolongation (Point point, Point p)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;

  struct { int x, y; } dp = {p.i - point.i, p.j - point.j};



   {
    if (dp.x == 0 && ((!is_leaf (neighbor(-1,0,0)) && neighbor(-1,0,0).neighbors && neighbor(-1,0,0).pid >= 0) || (!is_leaf (neighbor(1,0,0)) && neighbor(1,0,0).neighbors && neighbor(1,0,0).pid >= 0)))
      return true;
    if ((!is_leaf (neighbor(dp.x,0,0)) && neighbor(dp.x,0,0).neighbors && neighbor(dp.x,0,0).pid >= 0))
      return true;
  } 
#line 713
{
    if (dp.y == 0 && ((!is_leaf (neighbor(0,-1,0)) && neighbor(0,-1,0).neighbors && neighbor(0,-1,0).pid >= 0) || (!is_leaf (neighbor(0,1,0)) && neighbor(0,1,0).neighbors && neighbor(0,1,0).pid >= 0)))
      return true;
    if ((!is_leaf (neighbor(0,dp.y,0)) && neighbor(0,dp.y,0).neighbors && neighbor(0,dp.y,0).pid >= 0))
      return true;
  }
  return false;
}



static void append_pid (Array * pids, int pid)
{
  for (int i = 0, * p = (int *) pids->p; i < pids->len/sizeof(int); i++, p++)
    if (*p == pid)
      return;
  array_append (pids, &pid, sizeof(int));
}

static int locals_pids (Point point, Array * pids)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  if (is_leaf(cell)) {
    if (is_local(cell)) {
      Point p = point;
      {foreach_neighbor(1) {
 if ((cell.pid >= 0 && cell.pid != pid()) &&
     ((!is_leaf (cell) && cell.neighbors && cell.pid >= 0) || is_local_prolongation (point, p)))
   append_pid (pids, cell.pid);
 if ((!is_leaf (cell) && cell.neighbors && cell.pid >= 0))
   {foreach_child()
     if ((cell.pid >= 0 && cell.pid != pid()))
       append_pid (pids, cell.pid);end_foreach_child()}
      }end_foreach_neighbor()}
    }
  }
  else
    {foreach_neighbor(1) {
      if ((cell.pid >= 0 && cell.pid != pid()))
 append_pid (pids, cell.pid);
      if ((!is_leaf (cell) && cell.neighbors && cell.pid >= 0))
 {foreach_child()
   if ((cell.pid >= 0 && cell.pid != pid()))
     append_pid (pids, cell.pid);end_foreach_child()}
    }end_foreach_neighbor()}
  return pids->len/sizeof(int);
}

static int root_pids (Point point, Array * pids)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  {foreach_child()
    if ((cell.pid >= 0 && cell.pid != pid()))
      append_pid (pids, cell.pid);end_foreach_child()}
  return pids->len/sizeof(int);
}







static void rcv_pid_row (RcvPid * m, int l, int * row)
{
  for (int i = 0; i < npe(); i++)
    row[i] = 0;
  for (int i = 0; i < m->npid; i++) {
    Rcv * rcv = &m->rcv[i];
    if (l <= rcv->depth && rcv->halo[l].n > 0)
      row[rcv->pid] = rcv->halo[l].n;
  }
}

void check_snd_rcv_matrix (SndRcv * sndrcv, const char * name)
{
  int maxlevel = depth();
  mpi_all_reduce (maxlevel, MPI_INT, MPI_MAX);
  int * row = ((int *) pmalloc ((npe())*sizeof(int),__func__,__FILE__,0));
  for (int l = 0; l <= maxlevel; l++) {
    int status = 0;
    if (pid() == 0) {


      int ** send = matrix_new (npe(), npe(), sizeof(int));
      int ** receive = matrix_new (npe(), npe(), sizeof(int));
      rcv_pid_row (sndrcv->snd, l, row);
      MPI_Gather (row, npe(), MPI_INT, &send[0][0], npe(), MPI_INT, 0,
    MPI_COMM_WORLD);
      rcv_pid_row (sndrcv->rcv, l, row);
      MPI_Gather (row, npe(), MPI_INT, &receive[0][0], npe(), MPI_INT, 0,
    MPI_COMM_WORLD);

      int * astatus = ((int *) pmalloc ((npe())*sizeof(int),__func__,__FILE__,0));
      for (int i = 0; i < npe(); i++)
 astatus[i] = 0;
      for (int i = 0; i < npe(); i++)
 for (int j = 0; j < npe(); j++)
   if (send[i][j] != receive[j][i]) {
     fprintf (ferr, "%s: %d sends    %d to   %d at level %d\n",
       name, i, send[i][j], j, l);
     fprintf (ferr, "%s: %d receives %d from %d at level %d\n",
       name, j, receive[j][i], i, l);
     fflush (ferr);
     for (int k = i - 2; k <= i + 2; k++)
       if (k >= 0 && k < npe())
  astatus[k] = 1;
     for (int k = j - 2; k <= j + 2; k++)
       if (k >= 0 && k < npe())
  astatus[k] = 1;
   }
      MPI_Scatter (astatus, 1, MPI_INT, &status, 1, MPI_INT, 0, MPI_COMM_WORLD);
      pfree (astatus,__func__,__FILE__,0);

      matrix_free (send);
      matrix_free (receive);
    }
    else {
      rcv_pid_row (sndrcv->snd, l, row);
      MPI_Gather (row, npe(), MPI_INT, NULL, npe(), MPI_INT, 0, MPI_COMM_WORLD);
      rcv_pid_row (sndrcv->rcv, l, row);
      MPI_Gather (row, npe(), MPI_INT, NULL, npe(), MPI_INT, 0, MPI_COMM_WORLD);
      MPI_Scatter (NULL, 1, MPI_INT, &status, 1, MPI_INT, 0, MPI_COMM_WORLD);
    }
    if (status) {
      fprintf (ferr,
        "check_snd_rcv_matrix \"%s\" failed\n"
        "Calling debug_mpi(NULL)...\n"
        "Aborting...\n",
        name);
      fflush (ferr);
      debug_mpi (NULL);
      MPI_Abort (MPI_COMM_WORLD, -3);
    }
  }
  pfree (row,__func__,__FILE__,0);
}

static bool has_local_child (Point point)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  {foreach_child()
    if (is_local(cell))
      return true;end_foreach_child()}
  return false;
}

     
void mpi_boundary_update_buffers()
{tracing("mpi_boundary_update_buffers","/u/basilisk/src/grid/tree-mpi.h",0);
  if (npe() == 1)
    {end_tracing("mpi_boundary_update_buffers","/u/basilisk/src/grid/tree-mpi.h",0);return;}

  prof_start ("mpi_boundary_update_buffers");

  MpiBoundary * m = (MpiBoundary *) mpi_boundary;
  SndRcv * mpi_level = &m->mpi_level;
  SndRcv * mpi_level_root = &m->mpi_level_root;
  SndRcv * restriction = &m->restriction;

  snd_rcv_free (mpi_level);
  snd_rcv_free (mpi_level_root);
  snd_rcv_free (restriction);

  static const unsigned short used = 1 << user;
  {foreach_cell() {
    if (is_active(cell) && !is_border(cell))



      continue;

    if (cell.neighbors) {

      Array pids = {NULL, 0, 0};
      int n = locals_pids (point, &pids);
      if (n) {
 {foreach_child()
   if (is_local(cell))
     for (int i = 0, * p = (int *) pids.p; i < n; i++, p++)
       rcv_pid_append (mpi_level->snd, *p, point);end_foreach_child()}
 pfree (pids.p,__func__,__FILE__,0);
      }

      bool locals = false;
      if (is_leaf(cell)) {
 if ((cell.pid >= 0 && cell.pid != pid())) {
   Point p = point;
   {foreach_neighbor(1)
     if ((is_local(cell) &&
   ((!is_leaf (cell) && cell.neighbors && cell.pid >= 0) || is_local_prolongation (point, p))) ||
  is_root(point)) {
       locals = true; foreach_neighbor_break();
     }end_foreach_neighbor()}
 }
      }
      else
 {foreach_neighbor(1)
   if (is_local(cell) || is_root(point)) {
     locals = true; foreach_neighbor_break();
   }end_foreach_neighbor()}
      if (locals)
 {foreach_child()
   if ((cell.pid >= 0 && cell.pid != pid()))
            rcv_pid_append (mpi_level->rcv, cell.pid, point),
       cell.flags |= used;end_foreach_child()}


      if (!is_leaf(cell)) {

 if (is_local(cell)) {
   Array pids = {NULL, 0, 0};

   int n = root_pids (point, &pids);
   if (n) {
     {foreach_neighbor()
       for (int i = 0, * p = (int *) pids.p; i < n; i++, p++)
  if (cell.pid >= 0 && cell.pid != *p)
    rcv_pid_append (mpi_level_root->snd, *p, point);end_foreach_neighbor()}

     for (int i = 0, * p = (int *) pids.p; i < n; i++, p++)
       rcv_pid_append (restriction->snd, *p, point);
     pfree (pids.p,__func__,__FILE__,0);
   }
 }

 else if ((cell.pid >= 0 && cell.pid != pid())) {
   bool root = false;
   {foreach_child()
     if (is_local(cell)) {
       root = true; foreach_child_break();
     }end_foreach_child()}
   if (root) {
     int pid = cell.pid;
     {foreach_neighbor()
       if ((cell.pid >= 0 && cell.pid != pid()))
  rcv_pid_append (mpi_level_root->rcv, pid, point),
    cell.flags |= used;end_foreach_neighbor()}

     rcv_pid_append (restriction->rcv, pid, point);
   }
 }
      }
    }


    if (level > 0) {
      if (is_local(cell)) {

 Array pids = {NULL, 0, 0};
 if ((aparent(0,0,0).pid >= 0 && aparent(0,0,0).pid != pid()))
   append_pid (&pids, aparent(0,0,0).pid);
 int n = root_pids (parent, &pids);
 if (n) {
   for (int i = 0, * p = (int *) pids.p; i < n; i++, p++)
     rcv_pid_append (restriction->snd, *p, point);
   pfree (pids.p,__func__,__FILE__,0);
 }
      }
      else if ((cell.pid >= 0 && cell.pid != pid())) {

 if (is_local(aparent(0,0,0)) || has_local_child (parent))
   rcv_pid_append (restriction->rcv, cell.pid, point);
      }
    }
  }end_foreach_cell();}





  static const unsigned short keep = 1 << (user + 1);
  for (int l = depth(); l >= 0; l--)
    {foreach_cell()
      if (level == l) {
 if (level > 0 && (cell.pid < 0 || is_local(cell) || (cell.flags & used)))
   aparent(0,0,0).flags |= keep;
 if ((!is_leaf (cell) && cell.neighbors && cell.pid >= 0) && !(cell.flags & keep))
   coarsen_cell (point, NULL);
 cell.flags &= ~(used|keep);
 continue;
      }end_foreach_cell();}


  m->send->len = m->receive->len = 0;
  rcv_pid_append_pids (mpi_level->snd, m->send);
  rcv_pid_append_pids (mpi_level_root->snd, m->send);
  rcv_pid_append_pids (mpi_level->rcv, m->receive);
  rcv_pid_append_pids (mpi_level_root->rcv, m->receive);

  prof_stop();
#line 1015 "/u/basilisk/src/grid/tree-mpi.h"
end_tracing("mpi_boundary_update_buffers","/u/basilisk/src/grid/tree-mpi.h",0);}

     
void mpi_boundary_refine (scalar * list)
{tracing("mpi_boundary_refine","/u/basilisk/src/grid/tree-mpi.h",0);
  prof_start ("mpi_boundary_refine");

  MpiBoundary * mpi = (MpiBoundary *) mpi_boundary;


  Array * snd = mpi->send;
  MPI_Request r[2*snd->len/sizeof(int)];
  int nr = 0;
  for (int i = 0, * dest = snd->p; i < snd->len/sizeof(int); i++,dest++) {
    int len = ((Tree *)grid)->refined.n;
    MPI_Isend (&((Tree *)grid)->refined.n, 1, MPI_INT, *dest,
        (128), MPI_COMM_WORLD, &r[nr++]);
    if (len > 0)
      MPI_Isend (((Tree *)grid)->refined.p, sizeof(Index)/sizeof(int)*len,
   MPI_INT, *dest, (128), MPI_COMM_WORLD, &r[nr++]);
  }



  Array * rcv = mpi->receive;
  Cache rerefined = {NULL, 0, 0};
  for (int i = 0, * source = rcv->p; i < rcv->len/sizeof(int); i++,source++) {
    int len;
    mpi_recv_check (&len, 1, MPI_INT, *source, (128),
      MPI_COMM_WORLD, MPI_STATUS_IGNORE,
      "mpi_boundary_refine (len)");
    if (len > 0) {
      Index p[len];
      mpi_recv_check (p, sizeof(Index)/sizeof(int)*len,
        MPI_INT, *source, (128),
        MPI_COMM_WORLD, MPI_STATUS_IGNORE,
        "mpi_boundary_refine (p)");
      Cache refined = {p, len, len};
      {foreach_cache (refined)
 if (level <= depth() && allocated(0,0,0)) {
   if (is_leaf(cell)) {
     bool neighbors = false;
     {foreach_neighbor()
       if (allocated(0,0,0) && (is_active(cell) || is_local(aparent(0,0,0)))) {
  neighbors = true; foreach_neighbor_break();
       }end_foreach_neighbor()}

     if (neighbors)
       refine_cell (point, list, 0, &rerefined);
   }
 }end_foreach_cache();}
    }
  }


  if (nr)
    MPI_Waitall (nr, r, MPI_STATUSES_IGNORE);


  pfree (((Tree *)grid)->refined.p,__func__,__FILE__,0);
  ((Tree *)grid)->refined = rerefined;

  prof_stop();



  mpi_all_reduce (rerefined.n, MPI_INT, MPI_SUM);
  if (rerefined.n)
    mpi_boundary_refine (list);
  {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    _attribute[s.i].dirty = true;}}
end_tracing("mpi_boundary_refine","/u/basilisk/src/grid/tree-mpi.h",0);}

static void check_depth()
{
#line 1121 "/u/basilisk/src/grid/tree-mpi.h"
}

typedef struct {
  int refined, leaf;
} Remote;



     
void mpi_boundary_coarsen (int l, int too_fine)
{tracing("mpi_boundary_coarsen","/u/basilisk/src/grid/tree-mpi.h",0);
  if (npe() == 1)
    {end_tracing("mpi_boundary_coarsen","/u/basilisk/src/grid/tree-mpi.h",0);return;}

  check_depth();

  if (!(sizeof(Remote) == sizeof(double))) qassert ("/u/basilisk/src/grid/tree-mpi.h", 0, "sizeof(Remote) == sizeof(double)");

  scalar  remote=new_scalar("remote");
  {foreach_cell() {
    if (level == l) {
      if (is_local(cell)) {
 ((Remote *)&val(remote,0,0,0))->refined = (!is_leaf (cell) && cell.neighbors && cell.pid >= 0);
 ((Remote *)&val(remote,0,0,0))->leaf = is_leaf(cell);
      }
      else {
 ((Remote *)&val(remote,0,0,0))->refined = true;
 ((Remote *)&val(remote,0,0,0))->leaf = false;
      }
      continue;
    }
    if (is_leaf(cell))
      continue;
  }end_foreach_cell();}
  mpi_boundary_level (mpi_boundary, ((scalar[]){remote,{-1}}), l);

  {foreach_cell() {
    if (level == l) {
      if (!is_local(cell)) {
 if ((!is_leaf (cell) && cell.neighbors && cell.pid >= 0) && !((Remote *)&val(remote,0,0,0))->refined)
   coarsen_cell_recursive (point, NULL);
 else if (is_leaf(cell) && cell.neighbors && ((Remote *)&val(remote,0,0,0))->leaf) {
   int pid = cell.pid;
   {foreach_child()
     cell.pid = pid;end_foreach_child()}
 }
      }
      continue;
    }
    if (is_leaf(cell))
      continue;
  }end_foreach_cell();}

  check_depth();

  if (l > 0) {
    {foreach_cell() {
      if (level == l) {
 val(remote,0,0,0) = is_local(cell) ? cell.neighbors : 0;
 continue;
      }
      if (is_leaf(cell))
 continue;
    }end_foreach_cell();}
    mpi_boundary_level (mpi_boundary, ((scalar[]){remote,{-1}}), l);
    {foreach_cell() {
      if (level == l)
 if (!is_local(cell) && is_local(aparent(0,0,0)) && val(remote,0,0,0)) {
   aparent(0,0,0).flags &= ~too_fine;
   continue;
 }
      if (is_leaf(cell))
 continue;
    }end_foreach_cell();}
  }delete((scalar*)((scalar[]){remote,{-1}}));
end_tracing("mpi_boundary_coarsen","/u/basilisk/src/grid/tree-mpi.h",0);}

static void flag_border_cells()
{
  {foreach_cell() {
    if (is_active(cell)) {
      short flags = cell.flags & ~border;
      {foreach_neighbor() {
 if (!is_local(cell) || (level > 0 && !is_local(aparent(0,0,0)))) {
   flags |= border; foreach_neighbor_break();
 }

 if (is_refined_check())
   {foreach_child()
     if (!is_local(cell)) {
       flags |= border; foreach_child_break();
     }end_foreach_child()}
 if (flags & border)
   foreach_neighbor_break();
      }end_foreach_neighbor()}
      cell.flags = flags;
    }
    else {
      cell.flags &= ~border;

    }
    if (is_leaf(cell)) {
      if (cell.neighbors) {
 {foreach_child()
   cell.flags &= ~border;end_foreach_child()}
 if (is_border(cell)) {
   bool remote = false;
   {foreach_neighbor (2/2)
     if (!is_local(cell)) {
       remote = true; foreach_neighbor_break();
     }end_foreach_neighbor()}
   if (remote)
     {foreach_child()
       cell.flags |= border;end_foreach_child()}
 }
      }
      continue;
    }
  }end_foreach_cell();}
}

static int balanced_pid (long index, long nt, int nproc)
{
  long ne = max(1, nt/nproc), nr = nt % nproc;
  int pid = index < nr*(ne + 1) ?
    index/(ne + 1) :
    nr + (index - nr*(ne + 1))/ne;
  return min(nproc - 1, pid);
}


     
void mpi_partitioning()
{tracing("mpi_partitioning","/u/basilisk/src/grid/tree-mpi.h",0);
  prof_start ("mpi_partitioning");

  long nt = 0;
  
#if _OPENMP
  #undef OMP
  #define OMP(x)
#endif
{
#line 1258
foreach ()
    nt++;end_foreach();}
#if _OPENMP
  #undef OMP
  #define OMP(x) _Pragma(#x)
#endif



  
#line 1262
long i = 0;
  ((Tree *)grid)->dirty = true;
  {foreach_cell_post (is_active (cell))
    if (is_active (cell)) {
      if (is_leaf (cell)) {
 cell.pid = balanced_pid (i++, nt, npe());
 if (cell.neighbors > 0) {
   int pid = cell.pid;
   {foreach_child()
     cell.pid = pid;end_foreach_child()}
 }
 if (!is_local(cell))
   cell.flags &= ~active;
      }
      else {
 cell.pid = child(0,0,0).pid;
 bool inactive = true;
 {foreach_child()
   if (is_active(cell)) {
     inactive = false; foreach_child_break();
   }end_foreach_child()}
 if (inactive)
   cell.flags &= ~active;
      }
    }end_foreach_cell_post();}

  flag_border_cells();

  prof_stop();

  mpi_boundary_update_buffers();
end_tracing("mpi_partitioning","/u/basilisk/src/grid/tree-mpi.h",0);}

void restore_mpi (FILE * fp, scalar * list1)
{
  long index = 0, nt = 0, start = ftell (fp);
  scalar  size=new_scalar("size"), * list = list_concat (((scalar[]){size,{-1}}), list1);;
  long offset = sizeof(double)*list_len(list);


  static const unsigned short set = 1 << user;
  scalar * listm = is_constant(cm) ? NULL : (scalar *)((vector[]){fm,{{-1},{-1}}});
  {foreach_cell()
    if (balanced_pid (index, nt, npe()) <= pid()) {
      unsigned flags;
      if (fread (&flags, sizeof(unsigned), 1, fp) != 1) {
 fprintf (ferr, "restore(): error: expecting 'flags'\n");
 exit (1);
      }
      {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {
 double val;
 if (fread (&val, sizeof(double), 1, fp) != 1) {
   fprintf (ferr, "restore(): error: expecting scalar\n");
   exit (1);
 }
 if (s.i != INT_MAX)
   val(s,0,0,0) = val;
      }}}
      if (level == 0)
 nt = val(size,0,0,0);
      cell.pid = balanced_pid (index, nt, npe());
      cell.flags |= set;
      if (!(flags & leaf) && is_leaf(cell)) {
 if (balanced_pid (index + val(size,0,0,0) - 1, nt, npe()) < pid()) {
   fseek (fp, (sizeof(unsigned) + offset)*(val(size,0,0,0) - 1), SEEK_CUR);
   index += val(size,0,0,0);
   continue;
 }
 refine_cell (point, listm, 0, NULL);
      }
      index++;
      if (is_leaf(cell))
 continue;
    }end_foreach_cell();}


  fseek (fp, start, SEEK_SET);
  index = 0;
  {foreach_cell() {
    unsigned flags;
    if (fread (&flags, sizeof(unsigned), 1, fp) != 1) {
      fprintf (ferr, "restore(): error: expecting 'flags'\n");
      exit (1);
    }
    if (cell.flags & set)
      fseek (fp, offset, SEEK_CUR);
    else {
      {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {
 double val;
 if (fread (&val, sizeof(double), 1, fp) != 1) {
   fprintf (ferr, "restore(): error: expecting a scalar\n");
   exit (1);
 }
 if (s.i != INT_MAX)
   val(s,0,0,0) = val;
      }}}
      cell.pid = balanced_pid (index, nt, npe());
      if (is_leaf(cell) && cell.neighbors) {
 int pid = cell.pid;
 {foreach_child()
   cell.pid = pid;end_foreach_child()}
      }
    }
    if (!(flags & leaf) && is_leaf(cell)) {
      bool locals = false;
      {foreach_neighbor(1)
 if ((cell.flags & set) && (is_local(cell) || is_root(point))) {
   locals = true; foreach_neighbor_break();
 }end_foreach_neighbor()}
      if (locals)
 refine_cell (point, listm, 0, NULL);
      else {
 fseek (fp, (sizeof(unsigned) + offset)*(val(size,0,0,0) - 1), SEEK_CUR);
 index += val(size,0,0,0);
 continue;
      }
    }
    index++;
    if (is_leaf(cell))
      continue;
  }end_foreach_cell();}


  {foreach_cell_post (is_active (cell)) {
    cell.flags &= ~set;
    if (is_active (cell)) {
      if (is_leaf (cell)) {
 if (cell.neighbors > 0) {
   int pid = cell.pid;
   {foreach_child()
     cell.pid = pid;end_foreach_child()}
 }
 if (!is_local(cell))
   cell.flags &= ~active;
      }
      else if (!is_local(cell)) {
 bool inactive = true;
 {foreach_child()
   if (is_active(cell)) {
     inactive = false; foreach_child_break();
   }end_foreach_child()}
 if (inactive)
   cell.flags &= ~active;
      }
    }
  }end_foreach_cell_post();}

  flag_border_cells();

  mpi_boundary_update (list);
  pfree (list,__func__,__FILE__,0);delete((scalar*)((scalar[]){size,{-1}}));
}
#line 1435 "/u/basilisk/src/grid/tree-mpi.h"
     
double z_indexing (scalar index, bool leaves)
{tracing("z_indexing","/u/basilisk/src/grid/tree-mpi.h",0);



  scalar  size=new_scalar("size");
  subtree_size (size, leaves);






  double maxi = -1.;
  if (pid() == 0)
    {foreach_level(0)
      maxi = val(size,0,0,0) - 1.;end_foreach_level();}




  {foreach_level(0)
    val(index,0,0,0) = 0;end_foreach_level();}
  for (int l = 0; l < depth(); l++) {
    { Boundary ** _i = boundaries, * _b; while (_i && (_b = *_i++)) if (_b->restriction) _b->restriction (_b, ((scalar[]){index,{-1}}), l); };
    {foreach_cell() {
      if (level == l) {
 if (is_leaf(cell)) {
   if (is_local(cell) && cell.neighbors) {
     int i = val(index,0,0,0);
     {foreach_child()
       val(index,0,0,0) = i;end_foreach_child()}
   }
 }
 else {
   bool loc = is_local(cell);
   if (!loc)
     {foreach_child()
       if (is_local(cell)) {
  loc = true; foreach_child_break();
       }end_foreach_child()}
   if (loc) {
     int i = val(index,0,0,0) + !leaves;
     {foreach_child() {
       val(index,0,0,0) = i;
       i += val(size,0,0,0);
     }end_foreach_child()}
   }
 }
 continue;
      }
      if (is_leaf(cell))
 continue;
    }end_foreach_cell();}
  }
  { Boundary ** _i = boundaries, * _b; while (_i && (_b = *_i++)) if (_b->restriction) _b->restriction (_b, ((scalar[]){index,{-1}}), depth()); };

  {delete((scalar*)((scalar[]){size,{-1}}));{end_tracing("z_indexing","/u/basilisk/src/grid/tree-mpi.h",0);return maxi;}}delete((scalar*)((scalar[]){size,{-1}}));
end_tracing("z_indexing","/u/basilisk/src/grid/tree-mpi.h",0);}
#line 1687 "/u/basilisk/src/grid/tree.h"
#line 1 "grid/balance.h"
#line 1 "/u/basilisk/src/grid/balance.h"


typedef struct {
  short leaf, prolongation;
  int pid;
} NewPid;



#if TRASH
# define is_newpid() (!isnan(val(newpid,0,0,0)) && ((NewPid *)&val(newpid,0,0,0))->pid > 0)
#else
# define is_newpid() (((NewPid *)&val(newpid,0,0,0))->pid > 0)
#endif

Array * linear_tree (size_t size, scalar newpid)
{
  const unsigned short sent = 1 << user, next = 1 << (user + 1);
  Array * a = array_new();

  {foreach_cell_post_all (true)
    if (level > 0 && (cell.flags & (sent|next)))
      aparent(0,0,0).flags |= next;end_foreach_cell_post_all();}

  bool empty = true;
  {foreach_cell_all() {
    if (cell.flags & sent) {
      array_append (a, &cell, size);
      cell.flags &= ~sent;
      empty = false;
    }
    else {
      if (cell.pid >= 0 && ((NewPid *)&val(newpid,0,0,0))->leaf)
 if (!(is_leaf(cell))) qassert ("/u/basilisk/src/grid/balance.h", 0, "is_leaf(cell)");
      if (is_refined_check()) {


 bool prolo = false;
 {foreach_child()
   if (((NewPid *)&val(newpid,0,0,0))->prolongation)
     prolo = true;end_foreach_child()}
 if (prolo) {

   cell.flags |= leaf;
   array_append (a, &cell, sizeof(Cell));
   cell.flags &= ~leaf;
 }
 else
   array_append (a, &cell, sizeof(Cell));
      }
      else
 array_append (a, &cell, sizeof(Cell));
    }
    if (cell.flags & next)
      cell.flags &= ~next;
    else
      continue;
  }end_foreach_cell_all();}

  if (empty)
    a->len = 0;
  return a;
}

#define foreach_tree(t, size, list)\
{\
  const unsigned short _sent = 1 << user, _next = 1 << (user + 1);\
  scalar * _list = list;\
  char * _i = (char *) (t)->p;\
  foreach_cell_all() {\
    Cell * c = (Cell *) _i;\
    if (c->flags & _sent) {\
      _i += size;\

#line 74


#define end_foreach_tree()\
    }\
    else\
      _i += sizeof(Cell);\
    if (c->flags & _next) {\
      if (!(c->neighbors)) qassert ("/u/basilisk/src/grid/balance.h", 0, "c->neighbors");\
      if (!(c->flags & leaf) && is_leaf(cell) &&\
   (!is_newpid() || !((NewPid *)&val(newpid,0,0,0))->leaf))\
\
 refine_cell (point, _list, 0, NULL);\
      else if (!cell.neighbors)\
\
 alloc_children (point);\
    }\
    else\
      continue;\
  } end_foreach_cell_all();\
}\

#line 94


Array * neighborhood (scalar newpid, int nextpid, FILE * fp)
{
  const unsigned short sent = 1 << user;
  {foreach_cell() {

    bool root = false;
    if ((!is_local(cell) || ((NewPid *)&val(newpid,0,0,0))->pid - 1 != nextpid) && (!is_leaf (cell) && cell.neighbors && cell.pid >= 0)) {
      {foreach_child()
 if (is_local(cell) && ((NewPid *)&val(newpid,0,0,0))->pid - 1 == nextpid) {
   root = true; foreach_child_break();
 }end_foreach_child()}
      if (root && cell.pid != nextpid) {
 {foreach_neighbor()
   if (cell.pid != nextpid && is_newpid()) {
     if (fp)
       fprintf (fp, "%g %g %g %d %d root\n",
         x, y, z, ((NewPid *)&val(newpid,0,0,0))->pid - 1, cell.pid);
     cell.flags |= sent;
   }end_foreach_neighbor()}
      }
    }

    if ((is_local(cell) && ((NewPid *)&val(newpid,0,0,0))->pid - 1 == nextpid) || root) {
      {foreach_neighbor(1)
 if (cell.neighbors && cell.pid != nextpid)
   {foreach_child()
     if (cell.pid != nextpid && is_newpid()) {
       if (fp)
  fprintf (fp, "%g %g %g %d %d nextpid\n",
    x, y, z, ((NewPid *)&val(newpid,0,0,0))->pid - 1, cell.pid);
       cell.flags |= sent;
     }end_foreach_child()}end_foreach_neighbor()}
    }
    if (is_leaf(cell))
      continue;
  }end_foreach_cell();}

  return linear_tree (sizeof(Cell) + datasize, newpid);
}

static void send_tree (Array * a, int to, MPI_Request * r)
{
  MPI_Isend (&a->len, 1, MPI_LONG, to, (256), MPI_COMM_WORLD, &r[0]);
  if (a->len > 0) {
    MPI_Isend (a->p, a->len, MPI_BYTE, to, (256), MPI_COMM_WORLD, &r[1]);
    ((Tree *)grid)->dirty = true;
  }
}

static void receive_tree (int from, scalar newpid, FILE * fp)
{
  Array a;
  mpi_recv_check (&a.len, 1, MPI_LONG, from, (256),
    MPI_COMM_WORLD, MPI_STATUS_IGNORE, "receive_tree (len)");
  if (a.len > 0) {
    a.p = pmalloc (a.len,__func__,__FILE__,0);
    if (fp)
      fprintf (fp, "receiving %ld from %d\n", a.len, from);
    mpi_recv_check (a.p, a.len, MPI_BYTE, from, (256),
      MPI_COMM_WORLD, MPI_STATUS_IGNORE, "receive_tree (p)");

    {foreach_tree (&a, sizeof(Cell) + datasize, NULL) {
      memcpy (((char *)&cell) + sizeof(Cell), ((char *)c) + sizeof(Cell),
       datasize);
      if (!(((NewPid *)&val(newpid,0,0,0))->pid > 0)) qassert ("/u/basilisk/src/grid/balance.h", 0, "NEWPID()->pid > 0");
      if (fp)
 fprintf (fp, "%g %g %g %d %d %d %d %d %d recv\n",
   x, y, z, ((NewPid *)&val(newpid,0,0,0))->pid - 1, cell.pid,
   c->flags & leaf,
   cell.flags & leaf, from, ((NewPid *)&val(newpid,0,0,0))->leaf);
    }end_foreach_tree();}
    pfree (a.p,__func__,__FILE__,0);
    ((Tree *)grid)->dirty = true;
  }
}

static void wait_tree (Array * a, MPI_Request * r)
{
  MPI_Wait (&r[0], MPI_STATUS_IGNORE);
  if (a->len > 0)
    MPI_Wait (&r[1], MPI_STATUS_IGNORE);
}

static void check_flags()
{







}

struct {
  int min;
  bool leaves;

  int npe;
} mpi = {
  1,
  true
};

     
bool balance()
{tracing("balance","/u/basilisk/src/grid/balance.h",0);
  if (npe() == 1)
    {end_tracing("balance","/u/basilisk/src/grid/balance.h",0);return false;}

  if (!(sizeof(NewPid) == sizeof(double))) qassert ("/u/basilisk/src/grid/balance.h", 0, "sizeof(NewPid) == sizeof(double)");

  check_flags();

  long nl = 0, nt = 0;
  {foreach_cell() {
    if (is_local(cell)) {
      nt++;
      if (is_leaf(cell))
 nl++;
    }
    if (is_leaf(cell))
      continue;
  }end_foreach_cell();}

  grid->n = grid->tn = nl;
  grid->maxdepth = depth();
  long nmin = nl, nmax = nl;

  mpi_all_reduce (nmax, MPI_LONG, MPI_MAX);
  mpi_all_reduce (nmin, MPI_LONG, MPI_MIN);
  mpi_all_reduce (grid->tn, MPI_LONG, MPI_SUM);
  mpi_all_reduce (grid->maxdepth, MPI_INT, MPI_MAX);
  if (mpi.leaves)
    nt = grid->tn;
  else
    mpi_all_reduce (nt, MPI_LONG, MPI_SUM);

  long ne = max(1, nt/npe());

  if (ne < mpi.min) {
    mpi.npe = max(1, nt/mpi.min);
    ne = max(1, nt/mpi.npe);
  }
  else
    mpi.npe = npe();

  if (nmax - nmin <= 1)
    {end_tracing("balance","/u/basilisk/src/grid/balance.h",0);return false;}

  scalar  newpid=new_scalar("newpid");
  double zn = z_indexing (newpid, mpi.leaves);
  if (pid() == 0)
    if (!(zn + 1 == nt)) qassert ("/u/basilisk/src/grid/balance.h", 0, "zn + 1 == nt");

  FILE * fp = NULL;
#line 261 "/u/basilisk/src/grid/balance.h"
  bool next = false, prev = false;
  {foreach_cell_all() {
    if (is_local(cell)) {
      int pid = balanced_pid (val(newpid,0,0,0), nt, mpi.npe);
      pid = clamp (pid, cell.pid - 1, cell.pid + 1);
      if (pid == pid() + 1)
 next = true;
      else if (pid == pid() - 1)
 prev = true;
      ((NewPid *)&val(newpid,0,0,0))->pid = pid + 1;
      ((NewPid *)&val(newpid,0,0,0))->leaf = is_leaf(cell);
      ((NewPid *)&val(newpid,0,0,0))->prolongation = (!is_leaf(cell) && !cell.neighbors && cell.pid >= 0);
      if (fp)
 fprintf (fp, "%g %g %d %d newpid\n", x, y, ((NewPid *)&val(newpid,0,0,0))->pid - 1, cell.pid);
    }
    else
      val(newpid,0,0,0) = 0;
  }end_foreach_cell_all();}
  for (int l = 0; l <= depth(); l++)
    { Boundary ** _i = boundaries, * _b; while (_i && (_b = *_i++)) if (_b->level) _b->level (_b, ((scalar[]){newpid,{-1}}), l); };
#line 305 "/u/basilisk/src/grid/balance.h"
  Array * anext = next ? neighborhood (newpid, pid() + 1, fp) : array_new();
  Array * aprev = prev ? neighborhood (newpid, pid() - 1, fp) : array_new();

  if (fp)
    fflush (fp);

  check_flags();


  MPI_Request rprev[2], rnext[2];
  if (pid() > 0)
    send_tree (aprev, pid() - 1, rprev);
  if (pid() < npe() - 1)
    send_tree (anext, pid() + 1, rnext);


  if (pid() < npe() - 1)
    receive_tree (pid() + 1, newpid, fp);
  if (pid() > 0)
    receive_tree (pid() - 1, newpid, fp);


  if (pid() > 0)
    wait_tree (aprev, rprev);
  array_free (aprev);
  if (pid() < npe() - 1)
    wait_tree (anext, rnext);
  array_free (anext);

  if (fp)
    fflush (fp);


  int pid_changed = false;
  {foreach_cell_all() {
    if (cell.pid >= 0) {
      if (is_newpid()) {
 if (fp)
   fprintf (fp, "%g %g %g %d %d %d %d %d new\n",
     x, y, z, ((NewPid *)&val(newpid,0,0,0))->pid - 1, cell.pid,
     is_leaf(cell), cell.neighbors, ((NewPid *)&val(newpid,0,0,0))->leaf);
 if (cell.pid != ((NewPid *)&val(newpid,0,0,0))->pid - 1) {
   cell.pid = ((NewPid *)&val(newpid,0,0,0))->pid - 1;
   cell.flags &= ~(active|border);
   if (is_local(cell))
     cell.flags |= active;
   pid_changed = true;
 }
 if (((NewPid *)&val(newpid,0,0,0))->leaf && !is_leaf(cell) && cell.neighbors)
   coarsen_cell_recursive (point, NULL);
      }
      else if (level > 0 && ((NewPid *)&coarse(newpid,0,0,0))->leaf)
 cell.pid = aparent(0,0,0).pid;
    }

    if (!cell.neighbors && allocated_child(0,0,0)) {
      if (fp)
 fprintf (fp, "%g %g %g %d %d freechildren\n",
   x, y, z, ((NewPid *)&val(newpid,0,0,0))->pid - 1, cell.pid);
      free_children (point);
    }
  }end_foreach_cell_all();}

  if (((Tree *)grid)->dirty || pid_changed) {


    {foreach_cell_post (!is_leaf (cell))
      if (!is_leaf(cell) && !is_local(cell)) {
 unsigned short flags = cell.flags & ~active;
 {foreach_child()
   if (is_active(cell)) {
     flags |= active; foreach_child_break();
   }end_foreach_child()}
 cell.flags = flags;
      }end_foreach_cell_post();}

    flag_border_cells();
    pid_changed = true;
  }

  if (fp)
    fclose (fp);

  mpi_all_reduce (pid_changed, MPI_INT, MPI_MAX);
  if (pid_changed)
    mpi_boundary_update_buffers();

  {delete((scalar*)((scalar[]){newpid,{-1}}));{end_tracing("balance","/u/basilisk/src/grid/balance.h",0);return pid_changed;}}delete((scalar*)((scalar[]){newpid,{-1}}));
end_tracing("balance","/u/basilisk/src/grid/balance.h",0);}

void mpi_boundary_update (scalar * list)
{
  mpi_boundary_update_buffers();
  {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    _attribute[s.i].dirty = true;}}
  grid->tn = 0;
  boundary_internal ((scalar *)list, "/u/basilisk/src/grid/balance.h", 0);
  while (balance());
}
#line 1688 "/u/basilisk/src/grid/tree.h"
#else
void mpi_boundary_refine (scalar * list){}
void mpi_boundary_coarsen (int a, int b){}
void mpi_boundary_update (scalar * list) {
  {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    _attribute[s.i].dirty = true;}}
  boundary_internal ((scalar *)list, "/u/basilisk/src/grid/tree.h", 0);
}
#endif
#line 4 "/u/basilisk/src/grid/quadtree.h"

void quadtree_methods() {
  tree_methods();
}
#line 15 "cylinder-cpp.c"
#line 1 "cylinder.c"
#line 1 "/u/basilisk/src/ast/std/stdbool.h"
#include <stdbool.h>
#line 2 "cylinder.c"
#line 1 "/u/basilisk/src/ast/std/stdint.h"
#include <stdint.h>
#line 3 "cylinder.c"
#line 1 "embed.h"
#line 1 "/u/basilisk/src/embed.h"
#line 12 "/u/basilisk/src/embed.h"
#line 1 "fractions.h"
#line 1 "/u/basilisk/src/fractions.h"
#line 12 "/u/basilisk/src/fractions.h"
#line 1 "geometry.h"
#line 1 "/u/basilisk/src/geometry.h"
#line 28 "/u/basilisk/src/geometry.h"
double line_alpha (double c, coord n)
{
  double alpha, n1, n2;

  n1 = fabs (n.x); n2 = fabs (n.y);
  if (n1 > n2)
    do { double __tmp = n1; n1 = n2; n2 = __tmp; } while(0);

  c = clamp (c, 0., 1.);
  double v1 = n1/2.;
  if (c <= v1/n2)
    alpha = sqrt (2.*c*n1*n2);
  else if (c <= 1. - v1/n2)
    alpha = c*n2 + v1;
  else
    alpha = n1 + n2 - sqrt (2.*n1*n2*(1. - c));

  if (n.x < 0.)
    alpha += n.x;
  if (n.y < 0.)
    alpha += n.y;

  return alpha - (n.x + n.y)/2.;
}
#line 133 "/u/basilisk/src/geometry.h"
double line_area (double nx, double ny, double alpha)
{
  double a, v, area;

  alpha += (nx + ny)/2.;
  if (nx < 0.) {
    alpha -= nx;
    nx = - nx;
  }
  if (ny < 0.) {
    alpha -= ny;
    ny = - ny;
  }

  if (alpha <= 0.)
    return 0.;

  if (alpha >= nx + ny)
    return 1.;

  if (nx < 1e-10)
    area = alpha/ny;
  else if (ny < 1e-10)
    area = alpha/nx;
  else {
    v = sq(alpha);

    a = alpha - nx;
    if (a > 0.)
      v -= a*a;

    a = alpha - ny;
    if (a > 0.)
      v -= a*a;

    area = v/(2.*nx*ny);
  }

  return clamp (area, 0., 1.);
}
#line 237 "/u/basilisk/src/geometry.h"
double rectangle_fraction (coord n, double alpha, coord a, coord b)
{
  coord n1;
   {
    alpha -= n.x*(b.x + a.x)/2.;
    n1.x = n.x*(b.x - a.x);
  } 
#line 240
{
    alpha -= n.y*(b.y + a.y)/2.;
    n1.y = n.y*(b.y - a.y);
  }
  return line_area(n1.x, n1.y, alpha);
}
#line 262 "/u/basilisk/src/geometry.h"
int facets (coord n, double alpha, coord p[2])
{
  int i = 0;
  for (double s = -0.5; s <= 0.5; s += 1.)
    {
      if (fabs (n.y) > 1e-4 && i < 2) {
 double a = (alpha - s*n.x)/n.y;
 if (a >= -0.5 && a <= 0.5) {
   p[i].x = s;
   p[i++].y = a;
 }
      }
      
#line 267
if (fabs (n.x) > 1e-4 && i < 2) {
 double a = (alpha - s*n.y)/n.x;
 if (a >= -0.5 && a <= 0.5) {
   p[i].y = s;
   p[i++].x = a;
 }
      }}
  return i;
}
#line 352 "/u/basilisk/src/geometry.h"
double line_length_center (coord m, double alpha, coord * p)
{
  alpha += (m.x + m.y)/2.;

  coord n = m;
  
    if (n.x < 0.) {
      alpha -= n.x;
      n.x = - n.x;
    }
    
#line 358
if (n.y < 0.) {
      alpha -= n.y;
      n.y = - n.y;
    }

  p->x = p->y = p->z = 0.;

  if (alpha <= 0. || alpha >= n.x + n.y)
    return 0.;

  
    if (n.x < 1e-4) {
      p->x = 0.;
      p->y = (m.y < 0. ? 1. - alpha : alpha) - 0.5;
      return 1.;
    }
    
#line 369
if (n.y < 1e-4) {
      p->y = 0.;
      p->x = (m.x < 0. ? 1. - alpha : alpha) - 0.5;
      return 1.;
    }

  if (alpha >= n.x) {
    p->x += 1.;
    p->y += (alpha - n.x)/n.y;
  }
  else
    p->x += alpha/n.x;

  double ax = p->x, ay = p->y;
  if (alpha >= n.y) {
    p->y += 1.;
    ay -= 1.;
    p->x += (alpha - n.y)/n.x;
    ax -= (alpha - n.y)/n.x;
  }
  else {
    p->y += alpha/n.y;
    ay -= alpha/n.y;
  }

   {
    p->x /= 2.;
    p->x = clamp (p->x, 0., 1.);
    if (m.x < 0.)
      p->x = 1. - p->x;
    p->x -= 0.5;
  } 
#line 394
{
    p->y /= 2.;
    p->y = clamp (p->y, 0., 1.);
    if (m.y < 0.)
      p->y = 1. - p->y;
    p->y -= 0.5;
  }

  return sqrt (ax*ax + ay*ay);
}
#line 482 "/u/basilisk/src/geometry.h"
void line_center (coord m, double alpha, double a, coord * p)
{
  alpha += (m.x + m.y)/2.;

  coord n = m;
  
    if (n.x < 0.) {
      alpha -= n.x;
      n.x = - n.x;
    }
    
#line 488
if (n.y < 0.) {
      alpha -= n.y;
      n.y = - n.y;
    }

  p->z = 0.;
  if (alpha <= 0.) {
    p->x = p->y = -0.5;
    return;
  }

  if (alpha >= n.x + n.y) {
    p->x = p->y = 0.;
    return;
  }

  
    if (n.x < 1e-4) {
      p->x = 0.;
      p->y = sign(m.y)*(a/2. - 0.5);
      return;
    }
    
#line 505
if (n.y < 1e-4) {
      p->y = 0.;
      p->x = sign(m.x)*(a/2. - 0.5);
      return;
    }

  p->x = p->y = cube(alpha);

   {
    double b = alpha - n.x;
    if (b > 0.) {
      p->x -= sq(b)*(alpha + 2.*n.x);
      p->y -= cube(b);
    }
  } 
#line 513
{
    double b = alpha - n.y;
    if (b > 0.) {
      p->y -= sq(b)*(alpha + 2.*n.y);
      p->x -= cube(b);
    }
  }

   {
    p->x /= 6.*sq(n.x)*n.y*a;
    p->x = sign(m.x)*(p->x - 0.5);
  } 
#line 521
{
    p->y /= 6.*sq(n.y)*n.x*a;
    p->y = sign(m.y)*(p->y - 0.5);
  }
}
#line 13 "/u/basilisk/src/fractions.h"






#line 1 "myc2d.h"
#line 1 "/u/basilisk/src/myc2d.h"





coord mycs (Point point, scalar c)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  int ix;
  double c_t,c_b,c_r,c_l;
  double mx0,my0,mx1,my1,mm1,mm2;


  c_t = val(c,-1,1,0) + val(c,0,1,0) + val(c,1,1,0);
  c_b = val(c,-1,-1,0) + val(c,0,-1,0) + val(c,1,-1,0);
  c_r = val(c,1,-1,0) + val(c,1,0,0) + val(c,1,1,0);
  c_l = val(c,-1,-1,0) + val(c,-1,0,0) + val(c,-1,1,0);



  mx0 = 0.5*(c_l-c_r);
  my0 = 0.5*(c_b-c_t);


  if (fabs(mx0) <= fabs(my0)) {
    my0 = my0 > 0. ? 1. : -1.;
    ix = 1;
  }
  else {
    mx0 = mx0 > 0. ? 1. : -1.;
    ix = 0;
  }


  mm1 = val(c,-1,-1,0) + 2.0*val(c,-1,0,0) + val(c,-1,1,0);
  mm2 = val(c,1,-1,0) + 2.0*val(c,1,0,0) + val(c,1,1,0);
  mx1 = mm1 - mm2 + 1.e-30;
  mm1 = val(c,-1,-1,0) + 2.0*val(c,0,-1,0) + val(c,1,-1,0);
  mm2 = val(c,-1,1,0) + 2.0*val(c,0,1,0) + val(c,1,1,0);
  my1 = mm1 - mm2 + 1.e-30;


  if (ix) {
    mm1 = fabs(my1);
    mm1 = fabs(mx1)/mm1;
    if (mm1 > fabs(mx0)) {
      mx0 = mx1;
      my0 = my1;
    }
  }
  else {
    mm1 = fabs(mx1);
    mm1 = fabs(my1)/mm1;
    if (mm1 > fabs(my0)) {
      mx0 = mx1;
      my0 = my1;
    }
  }



  mm1 = fabs(mx0) + fabs(my0);
  coord n = {mx0/mm1, my0/mm1};

  return n;
}
#line 13 "/u/basilisk/src/fractions.h"






#line 1 "myc2d.h"
#line 1 "/u/basilisk/src/myc2d.h"





static void _stencil_mycs (Point point, scalar c)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;   
  
  
   


_stencil_val(c,-1,1,0); _stencil_val(c,0,1,0); _stencil_val(c,1,1,0); 


     
#line 14
_stencil_val(c,-1,-1,0); _stencil_val(c,0,-1,0); _stencil_val(c,1,-1,0); 
     _stencil_val(c,1,-1,0); _stencil_val(c,1,0,0); _stencil_val(c,1,1,0); 
     _stencil_val(c,-1,-1,0); _stencil_val(c,-1,0,0); _stencil_val(c,-1,1,0);
            
      
   
            
     
       
    



   
    





_stencil_val(c,-1,-1,0);_stencil_val(c,-1,0,0); _stencil_val(c,-1,1,0); 


     
  


      
#line 35
_stencil_val(c,1,-1,0);_stencil_val(c,1,0,0); _stencil_val(c,1,1,0);  
     
        _stencil_val(c,-1,-1,0);_stencil_val(c,0,-1,0); _stencil_val(c,1,-1,0);
       _stencil_val(c,-1,1,0);_stencil_val(c,0,1,0); _stencil_val(c,1,1,0);    
        
        
    
      
       
       
   
    
      
       


  
        
        
    
      
       
       
   



      
  

  
#line 64
return ;
}
#line 20 "/u/basilisk/src/fractions.h"
#line 41 "/u/basilisk/src/fractions.h"
void fraction_refine (Point point, scalar c)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;





  double cc = val(c,0,0,0);
  if (cc <= 0. || cc >= 1.)
    {foreach_child()
      val(c,0,0,0) = cc;end_foreach_child()}
  else {




    coord n = mycs (point, c);
    double alpha = line_alpha (cc, n);






    {foreach_child() {
      static const coord a = {0.,0.,0.}, b = {.5,.5,.5};
      coord nc;
      
 nc.x = child.x*n.x;
 
#line 69
nc.y = child.y*n.y;
      val(c,0,0,0) = rectangle_fraction (nc, alpha, a, b);
    }end_foreach_child()}
  }
}











static void alpha_refine (Point point, scalar alpha)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  vector n = _attribute[alpha.i].n;
  double alphac = 2.*val(alpha,0,0,0);
  coord m;
  
    m.x = val(n.x,0,0,0);
    
#line 91
m.y = val(n.y,0,0,0);
  {foreach_child() {
    val(alpha,0,0,0) = alphac;
    
      val(alpha,0,0,0) -= child.x*m.x/2.;
      
#line 95
val(alpha,0,0,0) -= child.y*m.y/2.;
  }end_foreach_child()}
}
#line 121 "/u/basilisk/src/fractions.h"
struct Fractions {
  scalar Phi;
  scalar c;
  vector s;
  double val;
};

     
void fractions (struct Fractions a)
{tracing("fractions","/u/basilisk/src/fractions.h",0);
  scalar Phi = a.Phi;
  scalar c = a.c;
  vector   s=(a.s).x.i?(a.s):new_face_vector("s");
  double val = a.val;
#line 145 "/u/basilisk/src/fractions.h"
  vector p;
  p.x = s.y; p.y = s.x;
#line 155 "/u/basilisk/src/fractions.h"
  foreach_face_stencil(){_stencil_is_face_y(){ {





_stencil_val(Phi,0,0,0);_stencil_val(Phi,1,0,0);{ {






      _stencil_val_a(p.x,0,0,0);_stencil_val(Phi,0,0,0);_stencil_val(Phi,0,0,0); _stencil_val(Phi,1,0,0);
_stencil_val(Phi,0,0,0);
 {_stencil_val_a(p.x,0,0,0); _stencil_val(p.x,0,0,0);   }     
         
    
#line 171
}
      








{_stencil_val_a(p.x,0,0,0);_stencil_val(Phi,0,0,0); _stencil_val(Phi,1,0,0);       }}





           
#line 180 "/u/basilisk/src/fractions.h"
    
  
}}end__stencil_is_face_y()
#line 155
_stencil_is_face_x(){ {





_stencil_val(Phi,0,0,0);_stencil_val(Phi,0,1,0);{ {






      _stencil_val_a(p.y,0,0,0);_stencil_val(Phi,0,0,0);_stencil_val(Phi,0,0,0); _stencil_val(Phi,0,1,0);
_stencil_val(Phi,0,0,0);
 {_stencil_val_a(p.y,0,0,0); _stencil_val(p.y,0,0,0);   }     
         
    
#line 171
}
      








{_stencil_val_a(p.y,0,0,0);_stencil_val(Phi,0,0,0); _stencil_val(Phi,0,1,0);       }}





           
#line 180 "/u/basilisk/src/fractions.h"
    
  
}}end__stencil_is_face_x()}end_foreach_face_stencil();
#line 155 "/u/basilisk/src/fractions.h"
  {foreach_face_generic(){is_face_y(){ {





    if ((val(Phi,0,0,0) - val)*(val(Phi,1,0,0) - val) < 0.) {






      val(p.x,0,0,0) = (val(Phi,0,0,0) - val)/(val(Phi,0,0,0) - val(Phi,1,0,0));
      if (val(Phi,0,0,0) < val)
 val(p.x,0,0,0) = 1. - val(p.x,0,0,0);
    }
#line 180 "/u/basilisk/src/fractions.h"
    else
      val(p.x,0,0,0) = (val(Phi,0,0,0) > val || val(Phi,1,0,0) > val);
  }}end_is_face_y()
#line 155
is_face_x(){ {





    if ((val(Phi,0,0,0) - val)*(val(Phi,0,1,0) - val) < 0.) {






      val(p.y,0,0,0) = (val(Phi,0,0,0) - val)/(val(Phi,0,0,0) - val(Phi,0,1,0));
      if (val(Phi,0,0,0) < val)
 val(p.y,0,0,0) = 1. - val(p.y,0,0,0);
    }
#line 180 "/u/basilisk/src/fractions.h"
    else
      val(p.y,0,0,0) = (val(Phi,0,0,0) > val || val(Phi,0,1,0) > val);
  }}end_is_face_x()}end_foreach_face_generic();}
#line 205 "/u/basilisk/src/fractions.h"
  scalar s_z = c;
  foreach_stencil()

  {    
#line 240 "/u/basilisk/src/fractions.h"
    
    
     { 
_stencil_val(p.y,0,0,0); _stencil_val(p.y,1,0,0);  
       
       
    
#line 245
} 
#line 242
{ 
_stencil_val(p.x,0,0,0); _stencil_val(p.x,0,1,0);  
       
       
    
#line 245
}





{
      {_stencil_val_a(s_z,0,0,0); _stencil_val(p.x,0,0,0); } 
{      





      
   






      
      for (int i = 0; i <= 1; i++)
 {
   {_stencil_val(p.x,0,i,0); _stencil_val(p.x,0,i,0); {       
     _stencil_val(p.x,0,i,0);_stencil_val(Phi,0,i,0); 
          
     
   }      }
   
#line 270
{_stencil_val(p.y,i,0,0); _stencil_val(p.y,i,0,0); {       
     _stencil_val(p.y,i,0,0);_stencil_val(Phi,i,0,0); 
          
     
   }      }}








{
 {_stencil_val_a(s_z,0,0,0);_stencil_val(p.x,0,0,0); _stencil_val(p.y,0,0,0);   }
{
 {_stencil_val_a(s_z,0,0,0);     } 
{



 _stencil_val_a(s_z,0,0,0);  

      }}}
#line 283 "/u/basilisk/src/fractions.h"
         
          
      
    







}}





       
    
  
#line 295
}end_foreach_stencil();
  {
#line 206
foreach()

  {
#line 240 "/u/basilisk/src/fractions.h"
    coord n;
    double nn = 0.;
     {
      n.x = val(p.y,0,0,0) - val(p.y,1,0,0);
      nn += fabs(n.x);
    } 
#line 242
{
      n.y = val(p.x,0,0,0) - val(p.x,0,1,0);
      nn += fabs(n.y);
    }





    if (nn == 0.)
      val(s_z,0,0,0) = val(p.x,0,0,0);
    else {





      
 n.x /= nn;
 
#line 260
n.y /= nn;






      double alpha = 0., ni = 0.;
      for (int i = 0; i <= 1; i++)
 {
   if (val(p.x,0,i,0) > 0. && val(p.x,0,i,0) < 1.) {
     double a = sign(val(Phi,0,i,0) - val)*(val(p.x,0,i,0) - 0.5);
     alpha += n.x*a + n.y*(i - 0.5);
     ni++;
   }
   
#line 270
if (val(p.y,i,0,0) > 0. && val(p.y,i,0,0) < 1.) {
     double a = sign(val(Phi,i,0,0) - val)*(val(p.y,i,0,0) - 0.5);
     alpha += n.y*a + n.x*(i - 0.5);
     ni++;
   }}
#line 283 "/u/basilisk/src/fractions.h"
      if (ni == 0)
 val(s_z,0,0,0) = max (val(p.x,0,0,0), val(p.y,0,0,0));
      else if (ni != 4)
 val(s_z,0,0,0) = line_area (n.x, n.y, alpha/ni);
      else {



 val(s_z,0,0,0) = 0.;

      }
    }
  }end_foreach();}{if(!(a.s).x.i)delete((scalar*)((vector[]){s,{{-1},{-1}}}));}
#line 347 "/u/basilisk/src/fractions.h"
end_tracing("fractions","/u/basilisk/src/fractions.h",0);}
#line 391 "/u/basilisk/src/fractions.h"
coord youngs_normal (Point point, scalar c)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  coord n;
  double nn = 0.;
  if (!(2 == 2)) qassert ("/u/basilisk/src/fractions.h", 0, "dimension == 2");
   {
    n.x = (val(c,-1,1,0) + 2.*val(c,-1,0,0) + val(c,-1,-1,0) -
    val(c,+1,1,0) - 2.*val(c,+1,0,0) - val(c,+1,-1,0));
    nn += fabs(n.x);
  } 
#line 396
{
    n.y = (val(c,1,-1,0) + 2.*val(c,0,-1,0) + val(c,-1,-1,0) -
    val(c,1,+1,0) - 2.*val(c,0,+1,0) - val(c,-1,+1,0));
    nn += fabs(n.y);
  }

  if (nn > 0.)
    {
      n.x /= nn;
      
#line 404
n.y /= nn;}
  else
    n.x = 1.;
  return n;
}





coord facet_normal (Point point, scalar c, vector s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  if (s.x.i >= 0) {
    coord n;
    double nn = 0.;
     {
      n.x = val(s.x,0,0,0) - val(s.x,1,0,0);
      nn += fabs(n.x);
    } 
#line 419
{
      n.y = val(s.y,0,0,0) - val(s.y,0,1,0);
      nn += fabs(n.y);
    }
    if (nn > 0.)
      {
 n.x /= nn;
 
#line 425
n.y /= nn;}
    else
      {
 n.x = 1./2;
 
#line 428
n.y = 1./2;}
    return n;
  }
  return mycs (point, c);
}






#line 414
static void _stencil_facet_normal (Point point, scalar c, vector s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  if (s.x.i >= 0) {    
    
    
     { 
_stencil_val(s.x,0,0,0); _stencil_val(s.x,1,0,0);  
       
       
    
#line 422
} 
#line 419
{ 
_stencil_val(s.y,0,0,0); _stencil_val(s.y,0,1,0);  
       
       
    
#line 422
}
      
   
       
    
       
  
    return ;
  } 
_stencil_mycs (point, c);
  
#line 431
return;
}
#line 441 "/u/basilisk/src/fractions.h"
     
void reconstruction (const scalar c, vector n, scalar alpha)
{tracing("reconstruction","/u/basilisk/src/fractions.h",0);
  foreach_stencil() {





_stencil_val(c,0,0,0); _stencil_val(c,0,0,0);{ {
      _stencil_val_a(alpha,0,0,0);  
      
 {_stencil_val_a(n.x,0,0,0);  }
 
#line 453
{_stencil_val_a(n.y,0,0,0);  }
    } 
{  






       _stencil_mycs (point, c);
      
 {_stencil_val_a(n.x,0,0,0);  }
 
#line 464
{_stencil_val_a(n.y,0,0,0);  }
      _stencil_val_a(alpha,0,0,0);_stencil_val(c,0,0,0);    
    }}





          
    
  
#line 467
}end_foreach_stencil();
  {
#line 444
foreach() {





    if (val(c,0,0,0) <= 0. || val(c,0,0,0) >= 1.) {
      val(alpha,0,0,0) = 0.;
      
 val(n.x,0,0,0) = 0.;
 
#line 453
val(n.y,0,0,0) = 0.;
    }
    else {






      coord m = mycs (point, c);
      
 val(n.x,0,0,0) = m.x;
 
#line 464
val(n.y,0,0,0) = m.y;
      val(alpha,0,0,0) = line_alpha (val(c,0,0,0), m);
    }
  }end_foreach();}
#line 476 "/u/basilisk/src/fractions.h"
  
    _attribute[n.x.i].refine = _attribute[n.x.i].prolongation = refine_injection;
    
#line 477
_attribute[n.y.i].refine = _attribute[n.y.i].prolongation = refine_injection;




  _attribute[alpha.i].n = n;
  _attribute[alpha.i].refine = _attribute[alpha.i].prolongation = alpha_refine;

end_tracing("reconstruction","/u/basilisk/src/fractions.h",0);}
#line 505 "/u/basilisk/src/fractions.h"
struct OutputFacets {
  scalar c;
  FILE * fp;
  vector s;
};

     
void output_facets (struct OutputFacets p)
{tracing("output_facets","/u/basilisk/src/fractions.h",0);
  scalar c = p.c;
  vector s = p.s;
  if (!p.fp) p.fp = fout;
  if (!s.x.i) s.x.i = -1;

  foreach_stencil()
    {_stencil_val(c,0,0,0); _stencil_val(c,0,0,0); {  
       _stencil_facet_normal (point, c, s);     
      _stencil_val(c,0,0,0); 

            
           
        
     
 
#line 538 "/u/basilisk/src/fractions.h"
    }        }end_foreach_stencil();

  {
#line 519
foreach()
    if (val(c,0,0,0) > 1e-6 && val(c,0,0,0) < 1. - 1e-6) {
      coord n = facet_normal (point, c, s);
      double alpha = line_alpha (val(c,0,0,0), n);

      coord segment[2];
      if (facets (n, alpha, segment) == 2)
 fprintf (p.fp, "%g %g\n%g %g\n\n",
   x + segment[0].x*Delta, y + segment[0].y*Delta,
   x + segment[1].x*Delta, y + segment[1].y*Delta);
#line 538 "/u/basilisk/src/fractions.h"
    }end_foreach();}

  fflush (p.fp);
end_tracing("output_facets","/u/basilisk/src/fractions.h",0);}







     
double interface_area (scalar c)
{tracing("interface_area","/u/basilisk/src/fractions.h",0);
  double area = 0.;
  foreach_stencil ()
    {_stencil_val(c,0,0,0); _stencil_val(c,0,0,0); {   
       _stencil_mycs (point, c);     
      _stencil_val(c,0,0,0); 
          
    }        }end_foreach_stencil();
  
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel reduction(+:area)){
#line 553
foreach ()
    if (val(c,0,0,0) > 1e-6 && val(c,0,0,0) < 1. - 1e-6) {
      coord n = mycs (point, c), p;
      double alpha = line_alpha (val(c,0,0,0), n);
      area += pow(Delta, 2 - 1)*line_length_center(n,alpha,&p);
    }end_foreach();mpi_all_reduce_array(&area,double,MPI_SUM,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
  
#line 559
{end_tracing("interface_area","/u/basilisk/src/fractions.h",0);return area;}
end_tracing("interface_area","/u/basilisk/src/fractions.h",0);}
#line 13 "/u/basilisk/src/embed.h"






scalar  cs={0};
vector  fs={{1},{2}};






#line 1 "embed-tree.h"
#line 1 "/u/basilisk/src/embed-tree.h"
#line 14 "/u/basilisk/src/embed-tree.h"
static void embed_fraction_refine (Point point, scalar cs)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  double cc = val(cs,0,0,0);





  if (cc <= 0. || cc >= 1.) {
    {foreach_child()
      val(cs,0,0,0) = cc;end_foreach_child()}
  }
  else {






    coord n = facet_normal (point, cs, fs);
    double alpha = line_alpha (cc, n);

    {foreach_child() {
      static const coord a = {0.,0.,0.}, b = {.5,.5,.5};
      coord nc;
      
 nc.x = child.x*n.x;
 
#line 40
nc.y = child.y*n.y;
      val(cs,0,0,0) = rectangle_fraction (nc, alpha, a, b);
    }end_foreach_child()}
  }
}








static void embed_face_fraction_refine_x (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  vector fs = _attribute[s.i].v;





  if (val(cs,0,0,0) <= 0. || val(cs,0,0,0) >= 1.) {





    for (int j = 0; j <= 1; j++)
      for (int k = 0; k <= 1; k++)
 fine(fs.x,1,j,k) = val(cs,0,0,0);
    for (int i = 0; i <= 1; i++)
      if (!(!is_leaf (neighbor(2*i-1,0,0)) && neighbor(2*i-1,0,0).neighbors && neighbor(2*i-1,0,0).pid >= 0) && neighbor(2*i-1,0,0).neighbors &&
   (is_local(cell) || is_local(neighbor(2*i-1,0,0))))
 for (int j = 0; j <= 1; j++)
   for (int k = 0; k <= 1; k++)
     fine(fs.x,2*i,j,k) = val(fs.x,i,0,0);
  }
  else {






    coord n = facet_normal (point, cs, fs);
    double alpha = line_alpha (val(cs,0,0,0), n);
#line 100 "/u/basilisk/src/embed-tree.h"
    if (2.*fabs(alpha) < fabs(n.y)) {
      double yc = alpha/n.y;
      int i = yc > 0.;
      fine(fs.x,1,1 - i,0) = n.y < 0. ? 1. - i : i;
      fine(fs.x,1,i,0) = n.y < 0. ? i - 2.*yc : 1. - i + 2.*yc;
    }
    else
      fine(fs.x,1,0,0) = fine(fs.x,1,1,0) = alpha > 0.;
#line 132 "/u/basilisk/src/embed-tree.h"
    for (int i = 0; i <= 1; i++)
      if (neighbor(2*i-1,0,0).neighbors &&
   (is_local(cell) || is_local(neighbor(2*i-1,0,0)))) {
 if (!(!is_leaf (neighbor(2*i-1,0,0)) && neighbor(2*i-1,0,0).neighbors && neighbor(2*i-1,0,0).pid >= 0)) {
   if (val(fs.x,i,0,0) <= 0. || val(fs.x,i,0,0) >= 1.)
     for (int j = 0; j <= 1; j++)
       for (int k = 0; k <= 1; k++)
  fine(fs.x,2*i,j,k) = val(fs.x,i,0,0);
   else {






     double a = val(fs.y,0,1,0) <= 0. || val(fs.y,2*i-1,1,0) <= 0. ||
       val(fs.y,0,0,0) >= 1. || val(fs.y,2*i-1,0,0) >= 1.;
     if ((2.*a - 1)*(val(fs.x,i,0,0) - 0.5) > 0.) {
       fine(fs.x,2*i,0,0) = a;
       fine(fs.x,2*i,1,0) = 2.*val(fs.x,i,0,0) - a;
     }
     else {
       fine(fs.x,2*i,0,0) = 2.*val(fs.x,i,0,0) + a - 1.;
       fine(fs.x,2*i,1,0) = 1. - a;
     }
#line 174 "/u/basilisk/src/embed-tree.h"
   }
 }




 for (int j = 0; j <= 1; j++)



     if (fine(fs.x,2*i,j,k) && !fine(cs,i,j,k))
       fine(fs.x,2*i,j,k) = 0.;
      }
  }
}

#line 53
static void embed_face_fraction_refine_y (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  vector fs = _attribute[s.i].v;





  if (val(cs,0,0,0) <= 0. || val(cs,0,0,0) >= 1.) {





    for (int j = 0; j <= 1; j++)
      for (int k = 0; k <= 1; k++)
 fine(fs.y,j,1,k) = val(cs,0,0,0);
    for (int i = 0; i <= 1; i++)
      if (!(!is_leaf (neighbor(0,2*i-1,0)) && neighbor(0,2*i-1,0).neighbors && neighbor(0,2*i-1,0).pid >= 0) && neighbor(0,2*i-1,0).neighbors &&
   (is_local(cell) || is_local(neighbor(0,2*i-1,0))))
 for (int j = 0; j <= 1; j++)
   for (int k = 0; k <= 1; k++)
     fine(fs.y,j,2*i,k) = val(fs.y,0,i,0);
  }
  else {






    coord n = facet_normal (point, cs, fs);
    double alpha = line_alpha (val(cs,0,0,0), n);
#line 100 "/u/basilisk/src/embed-tree.h"
    if (2.*fabs(alpha) < fabs(n.x)) {
      double yc = alpha/n.x;
      int i = yc > 0.;
      fine(fs.y,1 - i,1,0) = n.x < 0. ? 1. - i : i;
      fine(fs.y,i,1,0) = n.x < 0. ? i - 2.*yc : 1. - i + 2.*yc;
    }
    else
      fine(fs.y,0,1,0) = fine(fs.y,1,1,0) = alpha > 0.;
#line 132 "/u/basilisk/src/embed-tree.h"
    for (int i = 0; i <= 1; i++)
      if (neighbor(0,2*i-1,0).neighbors &&
   (is_local(cell) || is_local(neighbor(0,2*i-1,0)))) {
 if (!(!is_leaf (neighbor(0,2*i-1,0)) && neighbor(0,2*i-1,0).neighbors && neighbor(0,2*i-1,0).pid >= 0)) {
   if (val(fs.y,0,i,0) <= 0. || val(fs.y,0,i,0) >= 1.)
     for (int j = 0; j <= 1; j++)
       for (int k = 0; k <= 1; k++)
  fine(fs.y,j,2*i,k) = val(fs.y,0,i,0);
   else {






     double a = val(fs.x,1,0,0) <= 0. || val(fs.x,1,2*i-1,0) <= 0. ||
       val(fs.x,0,0,0) >= 1. || val(fs.x,0,2*i-1,0) >= 1.;
     if ((2.*a - 1)*(val(fs.y,0,i,0) - 0.5) > 0.) {
       fine(fs.y,0,2*i,0) = a;
       fine(fs.y,1,2*i,0) = 2.*val(fs.y,0,i,0) - a;
     }
     else {
       fine(fs.y,0,2*i,0) = 2.*val(fs.y,0,i,0) + a - 1.;
       fine(fs.y,1,2*i,0) = 1. - a;
     }
#line 174 "/u/basilisk/src/embed-tree.h"
   }
 }




 for (int j = 0; j <= 1; j++)



     if (fine(fs.y,j,2*i,k) && !fine(cs,j,i,k))
       fine(fs.y,j,2*i,k) = 0.;
      }
  }
}
#line 206 "/u/basilisk/src/embed-tree.h"




static inline void restriction_embed_linear (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;

  if (!val(cs,0,0,0)) {
    val(s,0,0,0) = 0.;
    return;
  }






  double val = 0., nv = 0.;
  for (int i = 0; i <= 1; i++)



      if (fine(cs,0,i,j) && fine(cs,1,!i,!j))
 val += (fine(s,0,i,j) + fine(s,1,!i,!j))/2., nv++;
  if (nv > 0.) {
    val(s,0,0,0) = val/nv;
    return;
  }





  coord p = {0.,0.,0.};
  {foreach_child()
    if (val(cs,0,0,0))
      p.x += x, p.y += y, p.z += z, val += val(s,0,0,0), nv++;end_foreach_child()}
  if (!(nv > 0.)) qassert ("/u/basilisk/src/embed-tree.h", 0, "nv > 0.");
  val(s,0,0,0) = val/nv;






  if (_attribute[s.i].embed_gradient && _attribute[s.i].boundary[0] != _attribute[s.i].boundary_homogeneous[0]) {
    coord o = {x,y,z}, g;
    _attribute[s.i].embed_gradient (point, s, &g);
    
      val(s,0,0,0) += (o.x - p.x/nv)*g.x;
      
#line 255
val(s,0,0,0) += (o.y - p.y/nv)*g.y;
  }
}
#line 268 "/u/basilisk/src/embed-tree.h"
static inline void refine_embed_linear (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  {foreach_child() {
    if (!val(cs,0,0,0))
      val(s,0,0,0) = 0.;
    else {
      if (!(coarse(cs,0,0,0))) qassert ("/u/basilisk/src/embed-tree.h", 0, "coarse(cs)");
      int i = (child.x + 1)/2, j = (child.y + 1)/2;

      if (coarse(fs.x,i,0,0) && coarse(fs.y,0,j,0) &&
   (coarse(cs,0,0,0) == 1. || coarse(cs,child.x,0,0) == 1. ||
    coarse(cs,0,child.y,0) == 1. || coarse(cs,child.x,child.y,0) == 1.)) {
 if (!(coarse(cs,child.x,0,0) && coarse(cs,0,child.y,0))) qassert ("/u/basilisk/src/embed-tree.h", 0, "coarse(cs,child.x) && coarse(cs,0,child.y)");
 if (coarse(fs.x,i,child.y,0) && coarse(fs.y,child.x,j,0)) {

   if (!(coarse(cs,child.x,child.y,0))) qassert ("/u/basilisk/src/embed-tree.h", 0, "coarse(cs,child.x,child.y)");
   val(s,0,0,0) = (9.*coarse(s,0,0,0) +
   3.*(coarse(s,child.x,0,0) + coarse(s,0,child.y,0)) +
   coarse(s,child.x,child.y,0))/16.;
 }
 else

   val(s,0,0,0) = (2.*coarse(s,0,0,0) + coarse(s,child.x,0,0) + coarse(s,0,child.y,0))/4.;
      }
      else if (coarse(cs,child.x,child.y,0) &&
        ((coarse(fs.x,i,0,0) && coarse(fs.y,child.x,j,0)) ||
  (coarse(fs.y,0,j,0) && coarse(fs.x,i,child.y,0)))) {

 val(s,0,0,0) = (3.*coarse(s,0,0,0) + coarse(s,child.x,child.y,0))/4.;
      }
#line 346 "/u/basilisk/src/embed-tree.h"
      else {

 val(s,0,0,0) = coarse(s,0,0,0);
  {
   if (coarse(fs.x,(child.x + 1)/2,0,0) && coarse(cs,child.x,0,0))
     val(s,0,0,0) += (coarse(s,child.x,0,0) - coarse(s,0,0,0))/4.;
   else if (coarse(fs.x,(- child.x + 1)/2,0,0) && coarse(cs,- child.x,0,0))
     val(s,0,0,0) -= (coarse(s,- child.x,0,0) - coarse(s,0,0,0))/4.;
 } 
#line 349
{
   if (coarse(fs.y,0,(child.y + 1)/2,0) && coarse(cs,0,child.y,0))
     val(s,0,0,0) += (coarse(s,0,child.y,0) - coarse(s,0,0,0))/4.;
   else if (coarse(fs.y,0,(- child.y + 1)/2,0) && coarse(cs,0,- child.y,0))
     val(s,0,0,0) -= (coarse(s,0,- child.y,0) - coarse(s,0,0,0))/4.;
 }
      }
    }
  }end_foreach_child()}
}
#line 369 "/u/basilisk/src/embed-tree.h"

void refine_embed_face_x (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  vector v = _attribute[s.i].v;
  for (int i = 0; i <= 1; i++)
    if (neighbor(2*i - 1,0,0).neighbors &&
 (is_local(cell) || is_local(neighbor(2*i - 1,0,0)))) {
      double g1 = val(fs.x,i,0,0) >= 1. && val(fs.x,i,+1,0) && val(fs.x,i,-1,0) ?
 (val(v.x,i,+1,0)/val(fs.x,i,+1,0) - val(v.x,i,-1,0)/val(fs.x,i,-1,0))/8. : 0.;
      double g2 = val(fs.x,i,0,0) >= 1. && val(fs.x,i,0,+1) && val(fs.x,i,0,-1) ?
 (val(v.x,i,0,+1)/val(fs.x,i,0,+1) - val(v.x,i,0,-1)/val(fs.x,i,0,-1))/8. : 0.;
      for (int j = 0; j <= 1; j++)
 for (int k = 0; k <= 1; k++)
   fine(v.x,2*i,j,k) = val(fs.x,i,0,0) ?
     fine(fs.x,2*i,j,k)*(val(v.x,i,0,0)/val(fs.x,i,0,0) +
    (2*j - 1)*g1 + (2*k - 1)*g2) : 0.;
    }
  if (is_local(cell)) {
    double g1 = (val(fs.x,0,+1,0) + val(fs.x,1,+1,0)) && (val(fs.x,0,-1,0) + val(fs.x,1,-1,0)) ?
      ((val(v.x,0,+1,0) + val(v.x,1,+1,0))/(val(fs.x,0,+1,0) + val(fs.x,1,+1,0)) -
       (val(v.x,0,-1,0) + val(v.x,1,-1,0))/(val(fs.x,0,-1,0) + val(fs.x,1,-1,0)))/8. : 0.;
    double g2 = (val(fs.x,1,0,+1) + val(fs.x,0,0,+1)) && (val(fs.x,1,0,-1) + val(fs.x,0,0,-1)) ?
      ((val(v.x,0,0,+1) + val(v.x,1,0,+1))/(val(fs.x,1,0,+1) + val(fs.x,0,0,+1)) -
       (val(v.x,0,0,-1) + val(v.x,1,0,-1))/(val(fs.x,1,0,-1) + val(fs.x,0,0,-1)))/8. : 0.;
    for (int j = 0; j <= 1; j++)
      for (int k = 0; k <= 1; k++)
 fine(v.x,1,j,k) = val(fs.x,0,0,0) + val(fs.x,1,0,0) ?
   fine(fs.x,1,j,k)*((val(v.x,0,0,0) + val(v.x,1,0,0))/(val(fs.x,0,0,0) + val(fs.x,1,0,0)) +
       (2*j - 1)*g1 + (2*k - 1)*g2) : 0.;
  }
}

#line 370
void refine_embed_face_y (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  vector v = _attribute[s.i].v;
  for (int i = 0; i <= 1; i++)
    if (neighbor(0,2*i - 1,0).neighbors &&
 (is_local(cell) || is_local(neighbor(0,2*i - 1,0)))) {
      double g1 = val(fs.y,0,i,0) >= 1. && val(fs.y,+1,i,0) && val(fs.y,-1,i,0) ?
 (val(v.y,+1,i,0)/val(fs.y,+1,i,0) - val(v.y,-1,i,0)/val(fs.y,-1,i,0))/8. : 0.;
      double g2 = val(fs.y,0,i,0) >= 1. && val(fs.y,0,i,+1) && val(fs.y,0,i,-1) ?
 (val(v.y,0,i,+1)/val(fs.y,0,i,+1) - val(v.y,0,i,-1)/val(fs.y,0,i,-1))/8. : 0.;
      for (int j = 0; j <= 1; j++)
 for (int k = 0; k <= 1; k++)
   fine(v.y,j,2*i,k) = val(fs.y,0,i,0) ?
     fine(fs.y,j,2*i,k)*(val(v.y,0,i,0)/val(fs.y,0,i,0) +
    (2*j - 1)*g1 + (2*k - 1)*g2) : 0.;
    }
  if (is_local(cell)) {
    double g1 = (val(fs.y,+1,0,0) + val(fs.y,+1,1,0)) && (val(fs.y,-1,0,0) + val(fs.y,-1,1,0)) ?
      ((val(v.y,+1,0,0) + val(v.y,+1,1,0))/(val(fs.y,+1,0,0) + val(fs.y,+1,1,0)) -
       (val(v.y,-1,0,0) + val(v.y,-1,1,0))/(val(fs.y,-1,0,0) + val(fs.y,-1,1,0)))/8. : 0.;
    double g2 = (val(fs.y,0,1,+1) + val(fs.y,0,0,+1)) && (val(fs.y,0,1,-1) + val(fs.y,0,0,-1)) ?
      ((val(v.y,0,0,+1) + val(v.y,0,1,+1))/(val(fs.y,0,1,+1) + val(fs.y,0,0,+1)) -
       (val(v.y,0,0,-1) + val(v.y,0,1,-1))/(val(fs.y,0,1,-1) + val(fs.y,0,0,-1)))/8. : 0.;
    for (int j = 0; j <= 1; j++)
      for (int k = 0; k <= 1; k++)
 fine(v.y,j,1,k) = val(fs.y,0,0,0) + val(fs.y,0,1,0) ?
   fine(fs.y,j,1,k)*((val(v.y,0,0,0) + val(v.y,0,1,0))/(val(fs.y,0,0,0) + val(fs.y,0,1,0)) +
       (2*j - 1)*g1 + (2*k - 1)*g2) : 0.;
  }
}
#line 28 "/u/basilisk/src/embed.h"
#line 68 "/u/basilisk/src/embed.h"

static inline double embed_face_gradient_x (Point point, scalar a, int i)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  int j = sign(val(fs.x,i,1,0) - val(fs.x,i,-1,0));
  if (!(val(cs,i,0,0) && val(cs,i-1,0,0))) qassert ("/u/basilisk/src/embed.h", 0, "cs[i] && cs[i-1]");
  if ((val(fs.x,i,j,0) > 0.5 && val(fs.y,i,j + (j < 0),0) && val(fs.y,i-1,j + (j < 0),0) && val(cs,i,j,0) && val(cs,i-1,j,0)))
    return ((1. + val(fs.x,i,0,0))*(val(a,i,0,0) - val(a,i-1,0,0)) +
     (1. - val(fs.x,i,0,0))*(val(a,i,j,0) - val(a,i-1,j,0)))/(2.*Delta);
  return (val(a,i,0,0) - val(a,i-1,0,0))/Delta;
}

#line 69
static inline double embed_face_gradient_y (Point point, scalar a, int i)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  int j = sign(val(fs.y,1,i,0) - val(fs.y,-1,i,0));
  if (!(val(cs,0,i,0) && val(cs,0,i-1,0))) qassert ("/u/basilisk/src/embed.h", 0, "cs[i] && cs[i-1]");
  if ((val(fs.y,j,i,0) > 0.5 && val(fs.x,j + (j < 0),i,0) && val(fs.x,j + (j < 0),i-1,0) && val(cs,j,i,0) && val(cs,j,i-1,0)))
    return ((1. + val(fs.y,0,i,0))*(val(a,0,i,0) - val(a,0,i-1,0)) +
     (1. - val(fs.y,0,i,0))*(val(a,j,i,0) - val(a,j,i-1,0)))/(2.*Delta);
  return (val(a,0,i,0) - val(a,0,i-1,0))/Delta;
}
#line 28 "/u/basilisk/src/embed.h"
#line 68 "/u/basilisk/src/embed.h"

static void _stencil_embed_face_gradient_x (Point point, scalar a, int i)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;    
   _stencil_val(fs.x,i,-1,0);_stencil_val(fs.x,i,1,0);
_stencil_val(cs,i,0,0); _stencil_val(cs,i-1,0,0);
_stencil_val(fs.x,i,o_stencil,0); _stencil_val(fs.y,i,o_stencil,0    ); _stencil_val(fs.y,i-1,o_stencil,0    ); _stencil_val(cs,i,o_stencil,0); _stencil_val(cs,i-1,o_stencil,0);
    { _stencil_val(fs.x,i,0,0);_stencil_val(a,i,0,0); _stencil_val(a,i-1,0,0); 
_stencil_val(fs.x,i,0,0);_stencil_val(a,i,o_stencil,0); _stencil_val(a,i-1,o_stencil,0);    
       
#line 75
}
_stencil_val(a,i,0,0); _stencil_val(a,i-1,0,0);  
      
         
  
#line 76
return  ;
}

#line 69
static void _stencil_embed_face_gradient_y (Point point, scalar a, int i)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;    
   _stencil_val(fs.y,-1,i,0);_stencil_val(fs.y,1,i,0);
_stencil_val(cs,0,i,0); _stencil_val(cs,0,i-1,0);
_stencil_val(fs.y,o_stencil,i,0); _stencil_val(fs.x,o_stencil,i,0    ); _stencil_val(fs.x,o_stencil,i-1,0    ); _stencil_val(cs,o_stencil,i,0); _stencil_val(cs,o_stencil,i-1,0);
    { _stencil_val(fs.y,0,i,0);_stencil_val(a,0,i,0); _stencil_val(a,0,i-1,0); 
_stencil_val(fs.y,0,i,0);_stencil_val(a,o_stencil,i,0); _stencil_val(a,o_stencil,i-1,0);    
       
#line 75
}
_stencil_val(a,0,i,0); _stencil_val(a,0,i-1,0);  
      
         
  
#line 76
return  ;
}


static inline double embed_face_value_x (Point point, scalar a, int i)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  int j = sign(val(fs.x,i,1,0) - val(fs.x,i,-1,0));
  return (val(fs.x,i,j,0) > 0.5 && val(fs.y,i,j + (j < 0),0) && val(fs.y,i-1,j + (j < 0),0) && val(cs,i,j,0) && val(cs,i-1,j,0)) ?
    ((1. + val(fs.x,i,0,0))*((val(a,i,0,0)*(1.5 + val(cs,i,0,0)) + val(a,i-1,0,0)*(1.5 + val(cs,i-1,0,0)))/ (val(cs,i,0,0) + val(cs,i-1,0,0) + 3.)) + (1. - val(fs.x,i,0,0))*((val(a,i,j,0)*(1.5 + val(cs,i,j,0)) + val(a,i-1,j,0)*(1.5 + val(cs,i-1,j,0)))/ (val(cs,i,j,0) + val(cs,i-1,j,0) + 3.)))/2. :
    ((val(a,i,0,0)*(1.5 + val(cs,i,0,0)) + val(a,i-1,0,0)*(1.5 + val(cs,i-1,0,0)))/ (val(cs,i,0,0) + val(cs,i-1,0,0) + 3.));
}

#line 80
static inline double embed_face_value_y (Point point, scalar a, int i)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  int j = sign(val(fs.y,1,i,0) - val(fs.y,-1,i,0));
  return (val(fs.y,j,i,0) > 0.5 && val(fs.x,j + (j < 0),i,0) && val(fs.x,j + (j < 0),i-1,0) && val(cs,j,i,0) && val(cs,j,i-1,0)) ?
    ((1. + val(fs.y,0,i,0))*((val(a,0,i,0)*(1.5 + val(cs,0,i,0)) + val(a,0,i-1,0)*(1.5 + val(cs,0,i-1,0)))/ (val(cs,0,i,0) + val(cs,0,i-1,0) + 3.)) + (1. - val(fs.y,0,i,0))*((val(a,j,i,0)*(1.5 + val(cs,j,i,0)) + val(a,j,i-1,0)*(1.5 + val(cs,j,i-1,0)))/ (val(cs,j,i,0) + val(cs,j,i-1,0) + 3.)))/2. :
    ((val(a,0,i,0)*(1.5 + val(cs,0,i,0)) + val(a,0,i-1,0)*(1.5 + val(cs,0,i-1,0)))/ (val(cs,0,i,0) + val(cs,0,i-1,0) + 3.));
}



#line 80
static void _stencil_embed_face_value_x (Point point, scalar a, int i)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;    
   _stencil_val(fs.x,i,-1,0);_stencil_val(fs.x,i,1,0);
_stencil_val(fs.x,i,o_stencil,0); _stencil_val(fs.y,i,o_stencil,0    ); _stencil_val(fs.y,i-1,o_stencil,0    ); _stencil_val(cs,i,o_stencil,0); _stencil_val(cs,i-1,o_stencil,0); 
_stencil_val(fs.x,i,0,0);_stencil_val(a,i,0,0); _stencil_val(cs,i,0,0); _stencil_val(a,i-1,0,0); _stencil_val(cs,i-1,0,0);_stencil_val(cs,i,0,0); _stencil_val(cs,i-1,0,0); _stencil_val(fs.x,i,0,0);_stencil_val(a,i,o_stencil,0); _stencil_val(cs,i,o_stencil,0); _stencil_val(a,i-1,o_stencil,0); _stencil_val(cs,i-1,o_stencil,0);_stencil_val(cs,i,o_stencil,0); _stencil_val(cs,i-1,o_stencil,0);
_stencil_val(a,i,0,0); _stencil_val(cs,i,0,0); _stencil_val(a,i-1,0,0); _stencil_val(cs,i-1,0,0);_stencil_val(cs,i,0,0); _stencil_val(cs,i-1,0,0);
  
#line 83
return                  
                    
    ;
}

#line 80
static void _stencil_embed_face_value_y (Point point, scalar a, int i)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;    
   _stencil_val(fs.y,-1,i,0);_stencil_val(fs.y,1,i,0);
_stencil_val(fs.y,o_stencil,i,0); _stencil_val(fs.x,o_stencil,i,0    ); _stencil_val(fs.x,o_stencil,i-1,0    ); _stencil_val(cs,o_stencil,i,0); _stencil_val(cs,o_stencil,i-1,0); 
_stencil_val(fs.y,0,i,0);_stencil_val(a,0,i,0); _stencil_val(cs,0,i,0); _stencil_val(a,0,i-1,0); _stencil_val(cs,0,i-1,0);_stencil_val(cs,0,i,0); _stencil_val(cs,0,i-1,0); _stencil_val(fs.y,0,i,0);_stencil_val(a,o_stencil,i,0); _stencil_val(cs,o_stencil,i,0); _stencil_val(a,o_stencil,i-1,0); _stencil_val(cs,o_stencil,i-1,0);_stencil_val(cs,o_stencil,i,0); _stencil_val(cs,o_stencil,i-1,0);
_stencil_val(a,0,i,0); _stencil_val(cs,0,i,0); _stencil_val(a,0,i-1,0); _stencil_val(cs,0,i-1,0);_stencil_val(cs,0,i,0); _stencil_val(cs,0,i-1,0);
  
#line 83
return                  
                    
    ;
}
#line 175 "/u/basilisk/src/embed.h"

#line 220 "/u/basilisk/src/embed.h"
static inline
double embed_geometry (Point point, coord * p, coord * n)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  *n = facet_normal (point, cs, fs);
  double alpha = line_alpha (val(cs,0,0,0), *n);
  double area = line_length_center(*n,alpha,p);
  normalize (n);
  return area;
}
#line 220 "/u/basilisk/src/embed.h"
static void 
_stencil_embed_geometry (Point point,_stencil_undefined  * p,_stencil_undefined  * n)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES; 
_stencil_facet_normal (point, cs, fs);     
   
  
#line 224
_stencil_val(cs,0,0,0);   
   
  
  return ;
}





static inline
double embed_area_center (Point point, double * x1, double * y1, double * z1)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  double area = 0.;
  if (val(cs,0,0,0) > 0. && val(cs,0,0,0) < 1.) {
    coord n, p;
    area = embed_geometry (point, &p, &n);
    *x1 += p.x*Delta, *y1 += p.y*Delta, *z1 += p.z*Delta;
  }
  return area;
}
#line 253 "/u/basilisk/src/embed.h"
double embed_interpolate (Point point, scalar s, coord p)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  if (!(2 == 2)) qassert ("/u/basilisk/src/embed.h", 0, "dimension == 2");
  int i = sign(p.x), j = sign(p.y);
  if (val(cs,i,0,0) && val(cs,0,j,0) && val(cs,i,j,0))

    return ((val(s,0,0,0)*(1. - fabs(p.x)) + val(s,i,0,0)*fabs(p.x))*(1. - fabs(p.y)) +
     (val(s,0,j,0)*(1. - fabs(p.x)) + val(s,i,j,0)*fabs(p.x))*fabs(p.y));
  else {


    double val = val(s,0,0,0);
     {
      int i = sign(p.x);
      if (val(cs,i,0,0))
 val += fabs(p.x)*(val(s,i,0,0) - val(s,0,0,0));
      else if (val(cs,-i,0,0))
 val += fabs(p.x)*(val(s,0,0,0) - val(s,-i,0,0));
    } 
#line 265
{
      int i = sign(p.y);
      if (val(cs,0,i,0))
 val += fabs(p.y)*(val(s,0,i,0) - val(s,0,0,0));
      else if (val(cs,0,-i,0))
 val += fabs(p.y)*(val(s,0,0,0) - val(s,0,-i,0));
    }
    return val;
  }
}
#line 253 "/u/basilisk/src/embed.h"
static void _stencil_embed_interpolate (Point point, scalar s,_stencil_undefined * p)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;         
      
  
_stencil_val(cs,o_stencil,0,0); _stencil_val(cs,0,o_stencil,0); _stencil_val(cs,o_stencil,o_stencil,0);{

    {_stencil_val(s,0,0,0);_stencil_val(s, o_stencil,0,0);
_stencil_val(s,0,o_stencil,0); _stencil_val(s,o_stencil,o_stencil,0);         
      
#line 260
} 
{  


     _stencil_val(s,0,0,0);
     {   
      
_stencil_val(cs,o_stencil,0,0);{
 {_stencil_val(s,o_stencil,0,0); _stencil_val(s,0,0,0);   } 
{_stencil_val(cs,o_stencil,0,0);
 {_stencil_val(s,0,0,0);_stencil_val(s, o_stencil,0,0);   } }}
       
      
    
#line 271
} 
#line 265
{   
      
_stencil_val(cs,0,o_stencil,0);{
 {_stencil_val(s,0,o_stencil,0); _stencil_val(s,0,0,0);   } 
{_stencil_val(cs,0,o_stencil,0);
 {_stencil_val(s,0,0,0);_stencil_val(s,0, o_stencil,0);   } }}
       
      
    
#line 271
} 
    
  }}
     
  

#line 274
}
#line 283 "/u/basilisk/src/embed.h"
struct Cleanup {
  scalar c;
  vector s;
  double smin;
  bool opposite;
};

     
int fractions_cleanup (struct Cleanup p)
{tracing("fractions_cleanup","/u/basilisk/src/embed.h",0);
  scalar c = p.c;
  vector s = p.s;







  int changed = 1, schanged = 0, i;
  for (i = 0; i < 100 && changed; i++) {




    foreach_face_stencil(){_stencil_is_face_x(){
      {_stencil_val(s.x,0,0,0);_stencil_val(c,0,0,0);_stencil_val(c,-1,0,0); _stencil_val(s.x,0,0,0);
 {_stencil_val_a(s.x,0,0,0);  }        }}end__stencil_is_face_x()
#line 308
_stencil_is_face_y(){
      {_stencil_val(s.y,0,0,0);_stencil_val(c,0,0,0);_stencil_val(c,0,-1,0); _stencil_val(s.y,0,0,0);
 {_stencil_val_a(s.y,0,0,0);  }        }}end__stencil_is_face_y()}end_foreach_face_stencil();




    {
#line 308
foreach_face_generic(){is_face_x(){
      if (val(s.x,0,0,0) && ((!val(c,0,0,0) || !val(c,-1,0,0)) || val(s.x,0,0,0) < p.smin))
 val(s.x,0,0,0) = 0.;}end_is_face_x()
#line 308
is_face_y(){
      if (val(s.y,0,0,0) && ((!val(c,0,0,0) || !val(c,0,-1,0)) || val(s.y,0,0,0) < p.smin))
 val(s.y,0,0,0) = 0.;}end_is_face_y()}end_foreach_face_generic();}

    changed = 0;
    foreach_stencil()
      {_stencil_val(c,0,0,0); _stencil_val(c,0,0,0); {   
 
  {
   for (int i = 0; i <= 1; i++)
     {_stencil_val(s.x,i,0,0);
          } 









_stencil_val(s.x,0,0,0); _stencil_val(s.x,1,0,0);
     {_stencil_val_a(c,0,0,0);   }
#line 329 "/u/basilisk/src/embed.h"
          
 
} 
#line 316
{
   for (int i = 0; i <= 1; i++)
     {_stencil_val(s.y,0,i,0);
          } 









_stencil_val(s.y,0,0,0); _stencil_val(s.y,0,1,0);
     {_stencil_val_a(c,0,0,0);   }
#line 329 "/u/basilisk/src/embed.h"
          
 
}
   







{_stencil_val_a(c,0,0,0);   }







    
      
#line 341
}      }end_foreach_stencil();
    
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel reduction(+:changed)){
#line 313
foreach()
      if (val(c,0,0,0) > 0. && val(c,0,0,0) < 1.) {
 int n = 0;
  {
   for (int i = 0; i <= 1; i++)
     if (val(s.x,i,0,0) > 0.)
       n++;
#line 329 "/u/basilisk/src/embed.h"
   if (p.opposite && val(s.x,0,0,0) == 0. && val(s.x,1,0,0) == 0.)
     val(c,0,0,0) = 0., changed++;
 } 
#line 316
{
   for (int i = 0; i <= 1; i++)
     if (val(s.y,0,i,0) > 0.)
       n++;
#line 329 "/u/basilisk/src/embed.h"
   if (p.opposite && val(s.y,0,0,0) == 0. && val(s.y,0,1,0) == 0.)
     val(c,0,0,0) = 0., changed++;
 }







 if (n < 2)
   val(c,0,0,0) = 0., changed++;
      }end_foreach();mpi_all_reduce_array(&changed,int,MPI_SUM,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}

    
#line 343
schanged += changed;
  }
  if (changed)
    fprintf (ferr, "WARNING: fractions_cleanup() did not converge after "
      "%d iterations\n", i);
  {end_tracing("fractions_cleanup","/u/basilisk/src/embed.h",0);return schanged;}
end_tracing("fractions_cleanup","/u/basilisk/src/embed.h",0);}
#line 373 "/u/basilisk/src/embed.h"

static inline double dirichlet_gradient_x (Point point, scalar s, scalar cs,
        coord n, coord p, double bc,
        double * coef)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  
    n.x = - n.x;
    
#line 379
n.y = - n.y;
  double d[2], v[2] = {HUGE,HUGE};
  bool defined = true;
  
    if (defined && !val(fs.x,(n.x > 0.),0,0))
      defined = false;
    
#line 383
if (defined && !val(fs.y,0,(n.y > 0.),0))
      defined = false;
  if (defined)
    for (int l = 0; l <= 1; l++) {
      int i = (l + 1)*sign(n.x);
      d[l] = (i - p.x)/n.x;
      double y1 = p.y + d[l]*n.y;
      int j = y1 > 0.5 ? 1 : y1 < -0.5 ? -1 : 0;
      y1 -= j;

      if (val(fs.x,i + (i < 0),j,0) && val(fs.y,i,j,0) && val(fs.y,i,j+1,0) &&
   val(cs,i,j-1,0) && val(cs,i,j,0) && val(cs,i,j+1,0))
 v[l] = ((((val(s,i,j-1,0)))*((y1) - 1.) + ((val(s,i,j+1,0)))*((y1) + 1.))*(y1)/2. - ((val(s,i,j,0)))*((y1) - 1.)*((y1) + 1.));
#line 417 "/u/basilisk/src/embed.h"
      else
 break;
    }
  if (v[0] == HUGE) {





    d[0] = max(1e-3, fabs(p.x/n.x));
    *coef = - 1./(d[0]*Delta);
    return bc/(d[0]*Delta);
  }





  *coef = 0.;
  if (v[1] != HUGE)
    return (d[1]*(bc - v[0])/d[0] - d[0]*(bc - v[1])/d[1])/((d[1] - d[0])*Delta);
  return (bc - v[0])/(d[0]*Delta);
}

#line 374
static inline double dirichlet_gradient_y (Point point, scalar s, scalar cs,
        coord n, coord p, double bc,
        double * coef)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  
    n.y = - n.y;
    
#line 379
n.x = - n.x;
  double d[2], v[2] = {HUGE,HUGE};
  bool defined = true;
  
    if (defined && !val(fs.y,0,(n.y > 0.),0))
      defined = false;
    
#line 383
if (defined && !val(fs.x,(n.x > 0.),0,0))
      defined = false;
  if (defined)
    for (int l = 0; l <= 1; l++) {
      int i = (l + 1)*sign(n.y);
      d[l] = (i - p.y)/n.y;
      double y1 = p.x + d[l]*n.x;
      int j = y1 > 0.5 ? 1 : y1 < -0.5 ? -1 : 0;
      y1 -= j;

      if (val(fs.y,j,i + (i < 0),0) && val(fs.x,j,i,0) && val(fs.x,j+1,i,0) &&
   val(cs,j-1,i,0) && val(cs,j,i,0) && val(cs,j+1,i,0))
 v[l] = ((((val(s,j-1,i,0)))*((y1) - 1.) + ((val(s,j+1,i,0)))*((y1) + 1.))*(y1)/2. - ((val(s,j,i,0)))*((y1) - 1.)*((y1) + 1.));
#line 417 "/u/basilisk/src/embed.h"
      else
 break;
    }
  if (v[0] == HUGE) {





    d[0] = max(1e-3, fabs(p.y/n.y));
    *coef = - 1./(d[0]*Delta);
    return bc/(d[0]*Delta);
  }





  *coef = 0.;
  if (v[1] != HUGE)
    return (d[1]*(bc - v[0])/d[0] - d[0]*(bc - v[1])/d[1])/((d[1] - d[0])*Delta);
  return (bc - v[0])/(d[0]*Delta);
}
#line 373 "/u/basilisk/src/embed.h"

static void _stencil_dirichlet_gradient_x (Point point, scalar s, scalar cs,
_stencil_undefined *
        
#line 375
n,_stencil_undefined * p,_stencil_undefined * bc,
_stencil_undefined 
        
#line 376
* coef)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;       
  
       
  
  
  
    {_stencil_val(fs.x,o_stencil,0,0  ); 
          }
    
#line 383
{_stencil_val(fs.y,0,o_stencil,0  ); 
          }
    
for (int l = 0; l <= 1; l++) {                         
       
         
      
      
        

_stencil_val(fs.x,    o_stencil,o_stencil,0); _stencil_val(fs.y,o_stencil,o_stencil,0); _stencil_val(fs.y,o_stencil,o_stencil,0);
   _stencil_val(cs,o_stencil,o_stencil,0); _stencil_val(cs,o_stencil,o_stencil,0); _stencil_val(cs,o_stencil,o_stencil,0);
#line 393
{
 
{_stencil_val(s,o_stencil,o_stencil,0);_stencil_val(s,o_stencil,o_stencil,0);_stencil_val(s,o_stencil,o_stencil,0);              }
 
#line 418
}

            
#line 417 "/u/basilisk/src/embed.h"
      
    
}         
     
   
     





   
     
  







return   ;
}

#line 374
static void _stencil_dirichlet_gradient_y (Point point, scalar s, scalar cs,
_stencil_undefined *
        
#line 375
n,_stencil_undefined * p,_stencil_undefined * bc,
_stencil_undefined 
        
#line 376
* coef)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;       
  
       
  
  
  
    {_stencil_val(fs.y,0,o_stencil,0  ); 
          }
    
#line 383
{_stencil_val(fs.x,o_stencil,0,0  ); 
          }
    
for (int l = 0; l <= 1; l++) {                         
       
         
      
      
        

_stencil_val(fs.y,o_stencil,    o_stencil,0); _stencil_val(fs.x,o_stencil,o_stencil,0); _stencil_val(fs.x,o_stencil,o_stencil,0);
   _stencil_val(cs,o_stencil,o_stencil,0); _stencil_val(cs,o_stencil,o_stencil,0); _stencil_val(cs,o_stencil,o_stencil,0);
#line 393
{
 
{_stencil_val(s,o_stencil,o_stencil,0);_stencil_val(s,o_stencil,o_stencil,0);_stencil_val(s,o_stencil,o_stencil,0);              }
 
#line 418
}

            
#line 417 "/u/basilisk/src/embed.h"
      
    
}         
     
   
     





   
     
  







return   ;
}

double dirichlet_gradient (Point point, scalar s, scalar cs,
      coord n, coord p, double bc, double * coef)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;

  
    if (fabs(n.x) >= fabs(n.y))
      return dirichlet_gradient_x (point, s, cs, n, p, bc, coef);
    
#line 446
if (fabs(n.y) >= fabs(n.x))
      return dirichlet_gradient_y (point, s, cs, n, p, bc, coef);
#line 457 "/u/basilisk/src/embed.h"
  return HUGE;
}


#line 441
static void _stencil_dirichlet_gradient (Point point, scalar s, scalar cs,
_stencil_undefined *
      
#line 442
n,_stencil_undefined * p,_stencil_undefined * bc,_stencil_undefined  * coef)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;

  
      
{ _stencil_dirichlet_gradient_x (point, s, cs,NULL ,NULL ,NULL ,NULL );}
      
#line 447
{ _stencil_dirichlet_gradient_y (point, s, cs,NULL ,NULL ,NULL ,NULL );}
       
#line 457 "/u/basilisk/src/embed.h"
  return ;
}

bid embed;
#line 469 "/u/basilisk/src/embed.h"
static inline
coord embed_gradient (Point point, vector u, coord p, coord n)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  coord dudn;
   {
    bool dirichlet;
    double vb = _attribute[u.x.i].boundary[embed] (point, point, u.x, &dirichlet);
    if (dirichlet) {
      double val;
      dudn.x = dirichlet_gradient (point, u.x, cs, n, p, vb, &val);
    }
    else
      dudn.x = vb;
    if (dudn.x == HUGE)
      dudn.x = 0.;
  } 
#line 473
{
    bool dirichlet;
    double vb = _attribute[u.y.i].boundary[embed] (point, point, u.y, &dirichlet);
    if (dirichlet) {
      double val;
      dudn.y = dirichlet_gradient (point, u.y, cs, n, p, vb, &val);
    }
    else
      dudn.y = vb;
    if (dudn.y == HUGE)
      dudn.y = 0.;
  }
  return dudn;
}
#line 469 "/u/basilisk/src/embed.h"
static void 
_stencil_embed_gradient (Point point, vector u,_stencil_undefined * p,_stencil_undefined * n)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES; 
  
   {   
    
    default_stencil ( point,((scalar[]){ u.x,{-1}}) ); 
{ 
       
_stencil_dirichlet_gradient (point, u.x, cs,NULL ,NULL ,NULL ,NULL ); 
      
    
#line 479
}   
     
    
        
     
       
  
#line 484
} 
#line 473
{   
    
    default_stencil ( point,((scalar[]){ u.y,{-1}}) ); 
{ 
       
_stencil_dirichlet_gradient (point, u.y, cs,NULL ,NULL ,NULL ,NULL ); 
      
    
#line 479
}   
     
    
        
     
       
  
#line 484
}
  return ;
}
#line 507 "/u/basilisk/src/embed.h"
     
void embed_force (scalar p, vector u, vector mu, coord * Fp, coord * Fmu)
{tracing("embed_force","/u/basilisk/src/embed.h",0);
  coord Fps = {0}, Fmus = {0};
  foreach_stencil ()
    {_stencil_val(cs,0,0,0); _stencil_val(cs,0,0,0); {    






      
       _stencil_embed_geometry (point,NULL ,NULL );   
            
      _stencil_embed_interpolate (point, p,NULL );
       
  
#line 533 "/u/basilisk/src/embed.h"
      if (constant(mu.x) != 0.) {      
 
  { 
_stencil_val(mu.x,0,0,0); _stencil_val(mu.x,1,0,0); 
     _stencil_val(fs.x,0,0,0); _stencil_val(fs.x,1,0,0); 
    
 
#line 538
} 
#line 535
{ 
_stencil_val(mu.y,0,0,0); _stencil_val(mu.y,0,1,0); 
     _stencil_val(fs.y,0,0,0); _stencil_val(fs.y,0,1,0); 
    
 
#line 538
}  
      
#line 598 "/u/basilisk/src/embed.h"
     
  _stencil_embed_gradient (point, u,NULL ,NULL );
  
        
      }
    }      }end_foreach_stencil();
  
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel  reduction(+:Fmus)reduction(+:Fps)){
#line 511
foreach ()
    if (val(cs,0,0,0) > 0. && val(cs,0,0,0) < 1.) {






      coord n, b;
      double area = embed_geometry (point, &b, &n);
      area *= pow (Delta, 2 - 1);
      double Fn = area*embed_interpolate (point, p, b);
      
 Fps.x += Fn*n.x;
 
#line 524
Fps.y += Fn*n.y;
#line 533 "/u/basilisk/src/embed.h"
      if (constant(mu.x) != 0.) {
 double mua = 0., fa = 0.;
  {
   mua += val(mu.x,0,0,0) + val(mu.x,1,0,0);
   fa += val(fs.x,0,0,0) + val(fs.x,1,0,0);
 } 
#line 535
{
   mua += val(mu.y,0,0,0) + val(mu.y,0,1,0);
   fa += val(fs.y,0,0,0) + val(fs.y,0,1,0);
 }
 mua /= fa;
#line 598 "/u/basilisk/src/embed.h"
 if (!(2 == 2)) qassert ("/u/basilisk/src/embed.h", 0, "dimension == 2");
 coord dudn = embed_gradient (point, u, b, n);
 
   Fmus.x -= area*mua*(dudn.x*(sq(n.x) + 1.) + dudn.y*n.x*n.y);
   
#line 601
Fmus.y -= area*mua*(dudn.y*(sq(n.y) + 1.) + dudn.x*n.y*n.x);
      }
    }end_foreach();mpi_all_reduce_array(&Fmus.x,double,MPI_SUM,2);mpi_all_reduce_array(&Fps.x,double,MPI_SUM,2);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}

  
#line 605
*Fp = Fps; *Fmu = Fmus;
end_tracing("embed_force","/u/basilisk/src/embed.h",0);}
#line 615 "/u/basilisk/src/embed.h"
double embed_vorticity (Point point, vector u, coord p, coord n)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;




  coord dudn = embed_gradient (point, u, p, n);
#line 632 "/u/basilisk/src/embed.h"
  return dudn.y*n.x - dudn.x*n.y;
}
#line 654 "/u/basilisk/src/embed.h"
double embed_flux (Point point, scalar s, vector mu, double * val)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;





  *val = 0.;
  if (val(cs,0,0,0) >= 1. || val(cs,0,0,0) <= 0.)
    return 0.;





  bool dirichlet;
  double grad = _attribute[s.i].boundary[embed] (point, point, s, &dirichlet);
  if (!grad && !dirichlet)
    return 0.;





  coord n = facet_normal (point, cs, fs), p;
  double alpha = line_alpha (val(cs,0,0,0), n);
  double area = line_length_center(n,alpha,&p);





  double coef = 0.;
  if (dirichlet) {
    normalize (&n);
    grad = dirichlet_gradient (point, s, cs, n, p, grad, &coef);
  }




  double mua = 0., fa = 0.;
   {
    mua += val(mu.x,0,0,0) + val(mu.x,1,0,0);
    fa += val(fs.x,0,0,0) + val(fs.x,1,0,0);
  } 
#line 696
{
    mua += val(mu.y,0,0,0) + val(mu.y,0,1,0);
    fa += val(fs.y,0,0,0) + val(fs.y,0,1,0);
  }
  *val = - mua/(fa + 1e-30)*grad*area/Delta;
  return - mua/(fa + 1e-30)*coef*area/Delta;;
}
#line 711 "/u/basilisk/src/embed.h"
#undef neumann
#define neumann(expr) (data ? embed_area_center (point, &x, &y, &z),\
        *((bool *)data) = false, (expr) :\
        Delta*(expr) + val(_s,0,0,0))\

#line 715

#undef neumann_homogeneous
#define neumann_homogeneous() (data ? *((bool *)data) = false, (0) :\
       val(_s,0,0,0))\

#line 719

#undef dirichlet
#define dirichlet(expr) (data ? embed_area_center (point, &x, &y, &z),\
        *((bool *)data) = true, (expr) :\
        2.*(expr) - val(_s,0,0,0))\

#line 724

#undef dirichlet_homogeneous
#define dirichlet_homogeneous() (data ? *((bool *)data) = true, (0) :\
         - val(_s,0,0,0))\

#line 728

#line 739 "/u/basilisk/src/embed.h"
static inline double bilinear_embed (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  if (!coarse(cs,0,0,0) || !coarse(cs,child.x,0,0))
    return coarse(s,0,0,0);

  if (!coarse(cs,0,child.y,0) || !coarse(cs,child.x,child.y,0))
    return coarse(s,0,0,0);







  return bilinear (point, s);
}
#line 786 "/u/basilisk/src/embed.h"
     
void update_tracer (scalar f, vector uf, vector flux, double dt)
{tracing("update_tracer","/u/basilisk/src/embed.h",0);
#line 798 "/u/basilisk/src/embed.h"
  scalar  e=new_scalar("e");
  foreach_stencil() {




_stencil_val(cs,0,0,0);{
      {_stencil_val_a(e,0,0,0);  } 






{_stencil_val(cs,0,0,0);{ {
      
 {_stencil_val_r(f,0,0,0);_stencil_val(flux.x,0,0,0); _stencil_val(flux.x,1,0,0);   }
 
#line 814
{_stencil_val_r(f,0,0,0);_stencil_val(flux.y,0,0,0); _stencil_val(flux.y,0,1,0);   }
      _stencil_val_a(e,0,0,0);  
    } 
#line 827
{   
      
      for (int i = 0; i <= 1; i++)
 {
   {_stencil_val(uf.x,i,0,0);
     {_stencil_val(uf.x,i,0,0);  }   }
   
#line 831
{_stencil_val(uf.y,0,i,0);
     {_stencil_val(uf.y,0,i,0);  }   }}     
      _stencil_val(cs,0,0,0);   




      
      
 { _stencil_val(flux.x,0,0,0); _stencil_val(flux.x,1,0,0);  }
 
#line 840
{ _stencil_val(flux.y,0,0,0); _stencil_val(flux.y,0,1,0);  }
_stencil_val(cs,0,0,0);






{ {
 _stencil_val_r(f,0,0,0);  
 _stencil_val_a(e,0,0,0);  
      } 







{
 _stencil_val_r(f,0,0,0);     
 
 {foreach_neighbor(1)
   {_stencil_val(cs,0,0,0);  }end_foreach_neighbor()}
 _stencil_val_a(e,0,0,0);_stencil_val(cs,0,0,0);    
      }}
        






         







      
    
#line 866
}}   
#line 827 "/u/basilisk/src/embed.h"
    
#line 866
}}




       






    
  
#line 867
}end_foreach_stencil();
  {
#line 799
foreach() {




    if (val(cs,0,0,0) <= 0.)
      val(e,0,0,0) = 0.;






    else if (val(cs,0,0,0) >= 1.) {
      
 val(f,0,0,0) += dt*(val(flux.x,0,0,0) - val(flux.x,1,0,0))/Delta;
 
#line 814
val(f,0,0,0) += dt*(val(flux.y,0,0,0) - val(flux.y,0,1,0))/Delta;
      val(e,0,0,0) = 0.;
    }
#line 827 "/u/basilisk/src/embed.h"
    else {
      double umax = 0.;
      for (int i = 0; i <= 1; i++)
 {
   if (fabs(val(uf.x,i,0,0)) > umax)
     umax = fabs(val(uf.x,i,0,0));
   
#line 831
if (fabs(val(uf.y,0,i,0)) > umax)
     umax = fabs(val(uf.y,0,i,0));}
      double dtmax = Delta*val(cs,0,0,0)/(umax + 1e-30);




      double F = 0.;
      
 F += val(flux.x,0,0,0) - val(flux.x,1,0,0);
 
#line 840
F += val(flux.y,0,0,0) - val(flux.y,0,1,0);
      F /= Delta*val(cs,0,0,0);






      if (dt <= dtmax) {
 val(f,0,0,0) += dt*F;
 val(e,0,0,0) = 0.;
      }







      else {
 val(f,0,0,0) += dtmax*F;
 double scs = 0.;
 {foreach_neighbor(1)
   scs += sq(val(cs,0,0,0));end_foreach_neighbor()}
 val(e,0,0,0) = (dt - dtmax)*F*val(cs,0,0,0)/scs;
      }
    }
  }end_foreach();}





  foreach_stencil() {   
    
    {foreach_neighbor(1)
      { _stencil_val(e,0,0,0); }end_foreach_neighbor()}
    _stencil_val_r(f,0,0,0); _stencil_val(cs,0,0,0); 
  }end_foreach_stencil();





  {
#line 873
foreach() {
    double se = 0.;
    {foreach_neighbor(1)
      se += val(e,0,0,0);end_foreach_neighbor()}
    val(f,0,0,0) += val(cs,0,0,0)*se;
  }end_foreach();}delete((scalar*)((scalar[]){e,{-1}}));
end_tracing("update_tracer","/u/basilisk/src/embed.h",0);}







static int metric_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i = 0);*ip=i;*tp=t;return ret;}      static int metric(const int i,const double t,Event *_ev){tracing("metric","/u/basilisk/src/embed.h",0);
{
  foreach_stencil()
    {_stencil_val_a(cs,0,0,0);  }end_foreach_stencil();
  {
#line 889
foreach()
    val(cs,0,0,0) = 1.;end_foreach();}
  foreach_face_stencil(){_stencil_is_face_x(){
    {_stencil_val_a(fs.x,0,0,0);  }}end__stencil_is_face_x()
#line 891
_stencil_is_face_y(){
    {_stencil_val_a(fs.y,0,0,0);  }}end__stencil_is_face_y()}end_foreach_face_stencil();
  {
#line 891
foreach_face_generic(){is_face_x(){
    val(fs.x,0,0,0) = 1.;}end_is_face_x()
#line 891
is_face_y(){
    val(fs.y,0,0,0) = 1.;}end_is_face_y()}end_foreach_face_generic();}

  _attribute[cs.i].refine = embed_fraction_refine;
#line 904 "/u/basilisk/src/embed.h"
  _attribute[cs.i].prolongation = fraction_refine;
  
    _attribute[fs.x.i].prolongation = embed_face_fraction_refine_x;
    
#line 906
_attribute[fs.y.i].prolongation = embed_face_fraction_refine_y;







  restriction (((scalar[]){cs, fs.x,fs.y,{-1}}));


  if (!(is_constant (cm) || cm.i == cs.i)) qassert ("/u/basilisk/src/embed.h", 0, "is_constant (cm) || cm.i == cs.i");

  cm = cs;
  fm = fs;
}{end_tracing("metric","/u/basilisk/src/embed.h",0);return 0;}end_tracing("metric","/u/basilisk/src/embed.h",0);}




static int defaults_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i = 0);*ip=i;*tp=t;return ret;}      static int defaults(const int i,const double t,Event *_ev){tracing("defaults","/u/basilisk/src/embed.h",0); {
  display ((struct _display){"draw_vof (c = 'cs', s = 'fs', filled = -1, "
    "fc = {0.5,0.5,0.5}, order = 2);"});
}{end_tracing("defaults","/u/basilisk/src/embed.h",0);return 0;}end_tracing("defaults","/u/basilisk/src/embed.h",0);}
#line 4 "cylinder.c"
#line 1 "navier-stokes/centered.h"
#line 1 "/u/basilisk/src/navier-stokes/centered.h"
#line 27 "/u/basilisk/src/navier-stokes/centered.h"
#line 1 "./run.h"
#line 1 "/u/basilisk/src/run.h"
#line 9 "/u/basilisk/src/run.h"
double dt = 1.;

#line 1 "./utils.h"
#line 1 "/u/basilisk/src/utils.h"







double DT = 1e10, CFL = 0.5;




struct {

  long nc;

  long tnc;

  double t;

  double speed;

  timer gt;
} perf;





void update_perf() {
  perf.nc += grid->n;
  perf.tnc += grid->tn;
  perf.t = timer_elapsed (perf.gt);
  perf.speed = perf.tnc/perf.t;
}






typedef struct {
  double cpu;
  double real;
  double speed;
  double min;
  double avg;
  double max;
  size_t tnc;
  long mem;
} timing;






timing timer_timing (timer t, int i, size_t tnc, double * mpi)
{
  timing s;
#if _MPI
  s.avg = mpi_time - t.tm;
#endif
  clock_t end = clock();
  s.cpu = ((double) (end - t.c))/CLOCKS_PER_SEC;
  s.real = timer_elapsed (t);
  if (tnc == 0) {
    double n = 0;
    
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel reduction(+:n)){
#line 69
foreach() n++;end_foreach();mpi_all_reduce_array(&n,double,MPI_SUM,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
    
#line 70
s.tnc = n;
    tnc = n*i;
  }
  else
    s.tnc = tnc;
#if _GNU_SOURCE
  struct rusage usage;
  getrusage (RUSAGE_SELF, &usage);
  s.mem = usage.ru_maxrss;
#else
  s.mem = 0;
#endif
#if _MPI
  if (mpi)
    MPI_Allgather (&s.avg, 1, MPI_DOUBLE, mpi, 1, MPI_DOUBLE, MPI_COMM_WORLD);
  s.max = s.min = s.avg;
  mpi_all_reduce (s.max, MPI_DOUBLE, MPI_MAX);
  mpi_all_reduce (s.min, MPI_DOUBLE, MPI_MIN);
  mpi_all_reduce (s.avg, MPI_DOUBLE, MPI_SUM);
  mpi_all_reduce (s.real, MPI_DOUBLE, MPI_SUM);
  mpi_all_reduce (s.mem, MPI_LONG, MPI_SUM);
  s.real /= npe();
  s.avg /= npe();
  s.mem /= npe();
#else
  s.min = s.max = s.avg = 0.;
#endif
  s.speed = s.real > 0. ? tnc/s.real : -1.;
  return s;
}




void timer_print (timer t, int i, size_t tnc)
{
  timing s = timer_timing (t, i, tnc, NULL);
  fprintf (fout,
    "\n# " "Quadtree"
    ", %d steps, %g CPU, %.4g real, %.3g points.step/s, %d var\n",
    i, s.cpu, s.real, s.speed, (int) (datasize/sizeof(double)));
#if _MPI
  fprintf (fout,
    "# %d procs, MPI: min %.2g (%.2g%%) "
    "avg %.2g (%.2g%%) max %.2g (%.2g%%)\n",
    npe(),
    s.min, 100.*s.min/s.real,
    s.avg, 100.*s.avg/s.real,
    s.max, 100.*s.max/s.real);
#endif
}







typedef struct {
  double avg, rms, max, volume;
} norm;

norm normf (scalar f)
{
  double avg = 0., rms = 0., max = 0., volume = 0.;
  foreach_stencil(
)
    {_stencil_val(f,0,0,0);_stencil_val(cm,0,0,0); {   
      _stencil_val(f,0,0,0);   
         
_stencil_val(cm,0,0,0); 
       _stencil_val(cm,0,0,0); 
       _stencil_val(cm,0,0,0); 
       
    
#line 143
}       }end_foreach_stencil();
  
#line 135
if(!is_constant(cm)){
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel  reduction(+:volume)
   reduction(+:rms) reduction(+:avg)reduction(max:max)){
#line 135
foreach(
)
    if (val(f,0,0,0) != HUGE && (sq(Delta)*val(cm,0,0,0)) > 0.) {
      double v = fabs(val(f,0,0,0));
      if (v > max) max = v;
      volume += (sq(Delta)*val(cm,0,0,0));
      avg += (sq(Delta)*val(cm,0,0,0))*v;
      rms += (sq(Delta)*val(cm,0,0,0))*sq(v);
    }end_foreach();mpi_all_reduce_array(&volume,double,MPI_SUM,1);mpi_all_reduce_array(&rms,double,MPI_SUM,1);mpi_all_reduce_array(&avg,double,MPI_SUM,1);mpi_all_reduce_array(&max,double,MPI_MAX,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 143
}else {double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);
  
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel  reduction(+:volume)
   reduction(+:rms) reduction(+:avg)reduction(max:max)){
#line 135
foreach(
)
    if (val(f,0,0,0) != HUGE && (sq(Delta)*_const_cm) > 0.) {
      double v = fabs(val(f,0,0,0));
      if (v > max) max = v;
      volume += (sq(Delta)*_const_cm);
      avg += (sq(Delta)*_const_cm)*v;
      rms += (sq(Delta)*_const_cm)*sq(v);
    }end_foreach();mpi_all_reduce_array(&volume,double,MPI_SUM,1);mpi_all_reduce_array(&rms,double,MPI_SUM,1);mpi_all_reduce_array(&avg,double,MPI_SUM,1);mpi_all_reduce_array(&max,double,MPI_MAX,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 143
}
  norm n;
  n.avg = volume ? avg/volume : 0.;
  n.rms = volume ? sqrt(rms/volume) : 0.;
  n.max = max;
  n.volume = volume;
  return n;
}





typedef struct {
  double min, max, sum, stddev, volume;
} stats;

stats statsf (scalar f)
{
  double min = 1e100, max = -1e100, sum = 0., sum2 = 0., volume = 0.;
  foreach_stencil(
)
    {_stencil_val(cm,0,0,0); _stencil_val(f,0,0,0); {
_stencil_val(cm,0,0,0); 
       _stencil_val(cm,0,0,0);_stencil_val(f,0,0,0); 
       _stencil_val(cm,0,0,0);_stencil_val(f,0,0,0); 
       _stencil_val(f,0,0,0); { _stencil_val(f,0,0,0); }
_stencil_val(f,0,0,0); { _stencil_val(f,0,0,0); }
         
         
    
#line 171
}      }end_foreach_stencil();
  
#line 163
if(!is_constant(cm)){
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel  reduction(min:min)
   reduction(max:max) reduction(+:volume) reduction(+:sum2)reduction(+:sum)){
#line 163
foreach(
)
    if ((sq(Delta)*val(cm,0,0,0)) > 0. && val(f,0,0,0) != HUGE) {
      volume += (sq(Delta)*val(cm,0,0,0));
      sum += (sq(Delta)*val(cm,0,0,0))*val(f,0,0,0);
      sum2 += (sq(Delta)*val(cm,0,0,0))*sq(val(f,0,0,0));
      if (val(f,0,0,0) > max) max = val(f,0,0,0);
      if (val(f,0,0,0) < min) min = val(f,0,0,0);
    }end_foreach();mpi_all_reduce_array(&min,double,MPI_MIN,1);mpi_all_reduce_array(&max,double,MPI_MAX,1);mpi_all_reduce_array(&volume,double,MPI_SUM,1);mpi_all_reduce_array(&sum2,double,MPI_SUM,1);mpi_all_reduce_array(&sum,double,MPI_SUM,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 171
}else {double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);
  
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel  reduction(min:min)
   reduction(max:max) reduction(+:volume) reduction(+:sum2)reduction(+:sum)){
#line 163
foreach(
)
    if ((sq(Delta)*_const_cm) > 0. && val(f,0,0,0) != HUGE) {
      volume += (sq(Delta)*_const_cm);
      sum += (sq(Delta)*_const_cm)*val(f,0,0,0);
      sum2 += (sq(Delta)*_const_cm)*sq(val(f,0,0,0));
      if (val(f,0,0,0) > max) max = val(f,0,0,0);
      if (val(f,0,0,0) < min) min = val(f,0,0,0);
    }end_foreach();mpi_all_reduce_array(&min,double,MPI_MIN,1);mpi_all_reduce_array(&max,double,MPI_MAX,1);mpi_all_reduce_array(&volume,double,MPI_SUM,1);mpi_all_reduce_array(&sum2,double,MPI_SUM,1);mpi_all_reduce_array(&sum,double,MPI_SUM,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 171
}
  stats s;
  s.min = min, s.max = max, s.sum = sum, s.volume = volume;
  if (volume > 0.)
    sum2 -= sum*sum/volume;
  s.stddev = sum2 > 0. ? sqrt(sum2/volume) : 0.;
  return s;
}
#line 187 "/u/basilisk/src/utils.h"
static double generic_limiter (double r, double beta)
{
  double v1 = min (r, beta), v2 = min (beta*r, 1.);
  v1 = max (0., v1);
  return max (v1, v2);
}

double minmod (double s0, double s1, double s2) {
  return s1 == s0 ? 0. : generic_limiter ((s2 - s1)/(s1 - s0), 1.)*(s1 - s0);
}

double superbee (double s0, double s1, double s2) {
  return s1 == s0 ? 0. : generic_limiter ((s2 - s1)/(s1 - s0), 2.)*(s1 - s0);
}

double sweby (double s0, double s1, double s2) {
  return s1 == s0 ? 0. : generic_limiter ((s2 - s1)/(s1 - s0), 1.5)*(s1 - s0);
}
#line 213 "/u/basilisk/src/utils.h"
double theta = 1.3;

double minmod2 (double s0, double s1, double s2)
{
  if (s0 < s1 && s1 < s2) {
    double d1 = theta*(s1 - s0), d2 = (s2 - s0)/2., d3 = theta*(s2 - s1);
    if (d2 < d1) d1 = d2;
    return min(d1, d3);
  }
  if (s0 > s1 && s1 > s2) {
    double d1 = theta*(s1 - s0), d2 = (s2 - s0)/2., d3 = theta*(s2 - s1);
    if (d2 > d1) d1 = d2;
    return max(d1, d3);
  }
  return 0.;
}
#line 237 "/u/basilisk/src/utils.h"
void gradients (scalar * f, vector * g)
{
  if (!(list_len(f) == vectors_len(g))) qassert ("/u/basilisk/src/utils.h", 0, "list_len(f) == vectors_len(g)");
  foreach_stencil() {
    scalar s; vector v;
    {vector*_i0=g;scalar*_i1= f;if(_i0)for(v=*_i0,s=*_i1;_i0->x.i>= 0;v=*++_i0,s=*++_i1){ {
      if (_attribute[s.i].gradient)
 { {

_stencil_val(fs.x,0,0,0);_stencil_val(fs.x,1,0,0);{
     {_stencil_val_a(v.x,0,0,0);  }

     
{_stencil_val_a(v.x,0,0,0);_stencil_val(s,-1,0,0); _stencil_val(s,0,0,0); _stencil_val(s,1,0,0);   }}

      
   
 
#line 251
} 
#line 244
{

_stencil_val(fs.y,0,0,0);_stencil_val(fs.y,0,1,0);{
     {_stencil_val_a(v.y,0,0,0);  }

     
{_stencil_val_a(v.y,0,0,0);_stencil_val(s,0,-1,0); _stencil_val(s,0,0,0); _stencil_val(s,0,1,0);   }}

      
   
 
#line 251
}}
      else
 { {

_stencil_val(fs.x,0,0,0);_stencil_val(fs.x,1,0,0);{
     {_stencil_val_a(v.x,0,0,0);  }

     
{_stencil_val_a(v.x,0,0,0);_stencil_val(s,1,0,0); _stencil_val(s,-1,0,0);   }}

      
   
 
#line 260
} 
#line 253
{

_stencil_val(fs.y,0,0,0);_stencil_val(fs.y,0,1,0);{
     {_stencil_val_a(v.y,0,0,0);  }

     
{_stencil_val_a(v.y,0,0,0);_stencil_val(s,0,1,0); _stencil_val(s,0,-1,0);   }}

      
   
 
#line 260
}}
    }}}
  }end_foreach_stencil();
  {
#line 240
foreach() {
    scalar s; vector v;
    {vector*_i0=g;scalar*_i1= f;if(_i0)for(v=*_i0,s=*_i1;_i0->x.i>= 0;v=*++_i0,s=*++_i1){ {
      if (_attribute[s.i].gradient)
 { {

   if (!val(fs.x,0,0,0) || !val(fs.x,1,0,0))
     val(v.x,0,0,0) = 0.;
   else

     val(v.x,0,0,0) = _attribute[s.i].gradient (val(s,-1,0,0), val(s,0,0,0), val(s,1,0,0))/Delta;
 } 
#line 244
{

   if (!val(fs.y,0,0,0) || !val(fs.y,0,1,0))
     val(v.y,0,0,0) = 0.;
   else

     val(v.y,0,0,0) = _attribute[s.i].gradient (val(s,0,-1,0), val(s,0,0,0), val(s,0,1,0))/Delta;
 }}
      else
 { {

   if (!val(fs.x,0,0,0) || !val(fs.x,1,0,0))
     val(v.x,0,0,0) = 0.;
   else

     val(v.x,0,0,0) = (val(s,1,0,0) - val(s,-1,0,0))/(2.*Delta);
 } 
#line 253
{

   if (!val(fs.y,0,0,0) || !val(fs.y,0,1,0))
     val(v.y,0,0,0) = 0.;
   else

     val(v.y,0,0,0) = (val(s,0,1,0) - val(s,0,-1,0))/(2.*Delta);
 }}
    }}}
  }end_foreach();}
}
#line 280 "/u/basilisk/src/utils.h"
void vorticity (const vector u, scalar omega)
{
  foreach_stencil()
    {_stencil_val_a(omega,0,0,0);_stencil_val(fm.x,1,0,0); _stencil_val(fm.x,0,0,0);_stencil_val(u.y,0,0,0);
        _stencil_val(fm.x,1,0,0);_stencil_val(u.y,1,0,0); _stencil_val(fm.x,0,0,0);_stencil_val(u.y,-1,0,0);
_stencil_val(fm.y,0,1,0); _stencil_val(fm.y,0,0,0);_stencil_val(u.x,0,0,0);
        _stencil_val(fm.y,0,0,0);_stencil_val(u.x,0,-1,0); _stencil_val(fm.y,0,1,0);_stencil_val(u.x,0,1,0);_stencil_val(cm,0,0,0);      
             
#line 286
}end_foreach_stencil();
  
#line 282
if(!is_constant(fm.x) && !is_constant(cm)){{foreach()
    val(omega,0,0,0) = ((val(fm.x,1,0,0) - val(fm.x,0,0,0))*val(u.y,0,0,0) +
        val(fm.x,1,0,0)*val(u.y,1,0,0) - val(fm.x,0,0,0)*val(u.y,-1,0,0) -
        (val(fm.y,0,1,0) - val(fm.y,0,0,0))*val(u.x,0,0,0) +
        val(fm.y,0,0,0)*val(u.x,0,-1,0) - val(fm.y,0,1,0)*val(u.x,0,1,0))/(2.*val(cm,0,0,0)*Delta + 1e-30);end_foreach();}}else if(is_constant(fm.x) && !is_constant(cm)){struct{double x,y;}_const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);
  {
#line 282
foreach()
    val(omega,0,0,0) = ((_const_fm.x - _const_fm.x)*val(u.y,0,0,0) +
        _const_fm.x*val(u.y,1,0,0) - _const_fm.x*val(u.y,-1,0,0) -
        (_const_fm.y - _const_fm.y)*val(u.x,0,0,0) +
        _const_fm.y*val(u.x,0,-1,0) - _const_fm.y*val(u.x,0,1,0))/(2.*val(cm,0,0,0)*Delta + 1e-30);end_foreach();}}else if(!is_constant(fm.x) && is_constant(cm)){double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);
  {
#line 282
foreach()
    val(omega,0,0,0) = ((val(fm.x,1,0,0) - val(fm.x,0,0,0))*val(u.y,0,0,0) +
        val(fm.x,1,0,0)*val(u.y,1,0,0) - val(fm.x,0,0,0)*val(u.y,-1,0,0) -
        (val(fm.y,0,1,0) - val(fm.y,0,0,0))*val(u.x,0,0,0) +
        val(fm.y,0,0,0)*val(u.x,0,-1,0) - val(fm.y,0,1,0)*val(u.x,0,1,0))/(2.*_const_cm*Delta + 1e-30);end_foreach();}}else {struct{double x,y;}_const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);
  {
#line 282
foreach()
    val(omega,0,0,0) = ((_const_fm.x - _const_fm.x)*val(u.y,0,0,0) +
        _const_fm.x*val(u.y,1,0,0) - _const_fm.x*val(u.y,-1,0,0) -
        (_const_fm.y - _const_fm.y)*val(u.x,0,0,0) +
        _const_fm.y*val(u.x,0,-1,0) - _const_fm.y*val(u.x,0,1,0))/(2.*_const_cm*Delta + 1e-30);end_foreach();}}
}





double change (scalar s, scalar sn)
{
  double max = 0.;
  foreach_stencil() {
_stencil_val(cm,0,0,0); {     
       _stencil_val(sn,0,0,0);_stencil_val(s,0,0,0);   
       
  
    }
       
    
#line 302
_stencil_val_a(sn,0,0,0); _stencil_val(s,0,0,0); 
  }end_foreach_stencil();
  
#line 296
if(!is_constant(cm)){
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel reduction(max:max)){
#line 296
foreach() {
    if ((sq(Delta)*val(cm,0,0,0)) > 0.) {
      double ds = fabs (val(s,0,0,0) - val(sn,0,0,0));
      if (ds > max)
 max = ds;
    }
    val(sn,0,0,0) = val(s,0,0,0);
  }end_foreach();mpi_all_reduce_array(&max,double,MPI_MAX,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 303
}else {double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);
  
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel reduction(max:max)){
#line 296
foreach() {
    if ((sq(Delta)*_const_cm) > 0.) {
      double ds = fabs (val(s,0,0,0) - val(sn,0,0,0));
      if (ds > max)
 max = ds;
    }
    val(sn,0,0,0) = val(s,0,0,0);
  }end_foreach();mpi_all_reduce_array(&max,double,MPI_MAX,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 303
}
  return max;
}





scalar lookup_field (const char * name)
{
  if (name)
    {scalar*_i=(scalar*)( all);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      if (!strcmp (_attribute[s.i].name, name))
 return s;}}
  return (scalar){-1};
}

vector lookup_vector (const char * name)
{
  if (name) {
    char component[strlen(name) + 3];
    strcpy (component, name);
    strcat (component, ".x");
    {scalar*_i=(scalar*)( all);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      if (!strcmp (_attribute[s.i].name, component))
 return _attribute[s.i].v;}}
  }
  return (vector){{-1}};
}
#line 340 "/u/basilisk/src/utils.h"
#define foreach_segment(_S,_p) {\
  coord t = {(_S)[1].x - (_S)[0].x, (_S)[1].y - (_S)[0].y};\
  double norm = sqrt(sq(t.x) + sq(t.y));\
  if (!(norm > 0.)) qassert ("/u/basilisk/src/utils.h", 0, "norm > 0.");\
  t.x = t.x/norm + 1e-6, t.y = t.y/norm - 1.5e-6;\
  double alpha = ((_S)[0].x*((_S)[1].y - (_S)[0].y) -\
    (_S)[0].y*((_S)[1].x - (_S)[0].x))/norm;\
  foreach()\
    if (fabs(t.y*x - t.x*y - alpha) < 0.708*Delta) {\
      coord _o = {x,y}, _p[2];\
      int _n = 0;\
 if (t.x)\
   for (int _i = -1; _i <= 1 && _n < 2; _i += 2) {\
     _p[_n].x = _o.x + _i*Delta/2.;\
     double a = (_p[_n].x - (_S)[0].x)/t.x;\
     _p[_n].y = (_S)[0].y + a*t.y;\
     if (fabs(_p[_n].y - _o.y) <= Delta/2.) {\
       a = clamp (a, 0., norm);\
       _p[_n].x = (_S)[0].x + a*t.x, _p[_n].y = (_S)[0].y + a*t.y;\
       if (fabs(_p[_n].x - _o.x) <= Delta/2. &&\
    fabs(_p[_n].y - _o.y) <= Delta/2.)\
  _n++;\
     }\
   }\
\
 if (t.y)\
   for (int _i = -1; _i <= 1 && _n < 2; _i += 2) {\
     _p[_n].y = _o.y + _i*Delta/2.;\
     double a = (_p[_n].y - (_S)[0].y)/t.y;\
     _p[_n].x = (_S)[0].x + a*t.x;\
     if (fabs(_p[_n].x - _o.x) <= Delta/2.) {\
       a = clamp (a, 0., norm);\
       _p[_n].y = (_S)[0].y + a*t.y, _p[_n].x = (_S)[0].x + a*t.x;\
       if (fabs(_p[_n].y - _o.y) <= Delta/2. &&\
    fabs(_p[_n].x - _o.x) <= Delta/2.)\
  _n++;\
     }\
   }\
\
      if (_n == 2) {\

#line 380

#line 410 "/u/basilisk/src/utils.h"
#define end_foreach_segment() } } end_foreach(); }




void fields_stats()
{
  fprintf (ferr, "# t = %g, fields = {", t);
  {scalar*_i=(scalar*)( all);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    fprintf (ferr, " %s", _attribute[s.i].name);}}
  fputs (" }\n", ferr);
  fprintf (ferr, "# %12s: %12s %12s %12s %12s\n",
    "name", "min", "avg", "stddev", "max");
  {scalar*_i=(scalar*)( all);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {
    stats ss = statsf (s);
    fprintf (ferr, "# %12s: %12g %12g %12g %12g\n",
      _attribute[s.i].name, ss.min, ss.sum/ss.volume, ss.stddev, ss.max);
  }}}
}

#line 1 "./output.h"
#line 1 "/u/basilisk/src/output.h"
#line 37 "/u/basilisk/src/output.h"
struct OutputField {
  scalar * list;
  FILE * fp;
  int n;
  bool linear;
  double box[2][2];
};

     
void output_field (struct OutputField p)
{tracing("output_field","/u/basilisk/src/output.h",0);
  if (!p.list) p.list = all;
  if (p.n == 0) p.n = N;
  if (!p.fp) p.fp = fout;
  p.n++;
  if (p.box[0][0] == 0. && p.box[0][1] == 0. &&
      p.box[1][0] == 0. && p.box[1][1] == 0.) {
    p.box[0][0] = X0; p.box[0][1] = Y0;
    p.box[1][0] = X0 + L0; p.box[1][1] = Y0 + L0;
  }

  boundary_internal ((scalar *)p.list, "/u/basilisk/src/output.h", 0);
  int len = list_len(p.list);
  double Delta = 0.999999*(p.box[1][0] - p.box[0][0])/(p.n - 1);
  int ny = (p.box[1][1] - p.box[0][1])/Delta + 1;
  double ** field = (double **) matrix_new (p.n, ny, len*sizeof(double));
  for (int i = 0; i < p.n; i++) {
    double x = Delta*i + p.box[0][0];
    for (int j = 0; j < ny; j++) {
      double y = Delta*j + p.box[0][1];
      if (p.linear) {
 int k = 0;
 {scalar*_i=(scalar*)( p.list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
   field[i][len*j + k++] = interpolate ((struct _interpolate){s, x, y});}}
      }
      else {
 Point point = locate ((struct _locate){x, y});int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
 int k = 0;
 {scalar*_i=(scalar*)( p.list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
   field[i][len*j + k++] = point.level >= 0 ? val(s,0,0,0) : HUGE;}}
      }
    }
  }

  if (pid() == 0) {
#if _MPI
    MPI_Reduce (MPI_IN_PLACE, field[0], len*p.n*ny, MPI_DOUBLE, MPI_MIN, 0,
  MPI_COMM_WORLD);
#endif
    fprintf (p.fp, "# 1:x 2:y");
    int i = 3;
    {scalar*_i=(scalar*)( p.list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      fprintf (p.fp, " %d:%s", i++, _attribute[s.i].name);}}
    fputc('\n', p.fp);
    for (int i = 0; i < p.n; i++) {
      double x = Delta*i + p.box[0][0];
      for (int j = 0; j < ny; j++) {
 double y = Delta*j + p.box[0][1];

 fprintf (p.fp, "%g %g", x, y);
 int k = 0;
 {scalar*_i=(scalar*)( p.list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
   fprintf (p.fp, " %g", field[i][len*j + k++]);}}
 fputc ('\n', p.fp);
      }
      fputc ('\n', p.fp);
    }
    fflush (p.fp);
  }
#if _MPI
  else
    MPI_Reduce (field[0], NULL, len*p.n*ny, MPI_DOUBLE, MPI_MIN, 0,
  MPI_COMM_WORLD);
#endif

  matrix_free (field);
end_tracing("output_field","/u/basilisk/src/output.h",0);}
#line 141 "/u/basilisk/src/output.h"
struct OutputMatrix {
  scalar f;
  FILE * fp;
  int n;
  bool linear;
};

     
void output_matrix (struct OutputMatrix p)
{tracing("output_matrix","/u/basilisk/src/output.h",0);
  if (p.n == 0) p.n = N;
  if (!p.fp) p.fp = fout;
  if (p.linear) {
    scalar f = p.f;
    boundary_internal ((scalar *)((scalar[]){f,{-1}}), "/u/basilisk/src/output.h", 0);
  }
  float fn = p.n;
  float Delta = (float) L0/fn;
  fwrite (&fn, sizeof(float), 1, p.fp);
  for (int j = 0; j < p.n; j++) {
    float yp = (float) (Delta*j + X0 + Delta/2.);
    fwrite (&yp, sizeof(float), 1, p.fp);
  }
  for (int i = 0; i < p.n; i++) {
    float xp = (float) (Delta*i + X0 + Delta/2.);
    fwrite (&xp, sizeof(float), 1, p.fp);
    for (int j = 0; j < p.n; j++) {
      float yp = (float)(Delta*j + Y0 + Delta/2.), v;
      if (p.linear)
 v = interpolate ((struct _interpolate){p.f, xp, yp});
      else {
 Point point = locate ((struct _locate){xp, yp});int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
 if (!(point.level >= 0)) qassert ("/u/basilisk/src/output.h", 0, "point.level >= 0");
 v = val(p.f,0,0,0);
      }
      fwrite (&v, sizeof(float), 1, p.fp);
    }
  }
  fflush (p.fp);
end_tracing("output_matrix","/u/basilisk/src/output.h",0);}
#line 189 "/u/basilisk/src/output.h"
typedef void (* colormap) (double cmap[127][3]);

void jet (double cmap[127][3])
{
  for (int i = 0; i < 127; i++) {
    cmap[i][0] =
      i <= 46 ? 0. :
      i >= 111 ? -0.03125*(i - 111) + 1. :
      i >= 78 ? 1. :
      0.03125*(i - 46);
    cmap[i][1] =
      i <= 14 || i >= 111 ? 0. :
      i >= 79 ? -0.03125*(i - 111) :
      i <= 46 ? 0.03125*(i - 14) :
      1.;
    cmap[i][2] =
      i >= 79 ? 0. :
      i >= 47 ? -0.03125*(i - 79) :
      i <= 14 ? 0.03125*(i - 14) + 1.:
      1.;
  }
}

void cool_warm (double cmap[127][3])
{






  static double basemap[33][3] = {
    {0.2298057, 0.298717966, 0.753683153},
    {0.26623388, 0.353094838, 0.801466763},
    {0.30386891, 0.406535296, 0.84495867},
    {0.342804478, 0.458757618, 0.883725899},
    {0.38301334, 0.50941904, 0.917387822},
    {0.424369608, 0.558148092, 0.945619588},
    {0.46666708, 0.604562568, 0.968154911},
    {0.509635204, 0.648280772, 0.98478814},
    {0.552953156, 0.688929332, 0.995375608},
    {0.596262162, 0.726149107, 0.999836203},
    {0.639176211, 0.759599947, 0.998151185},
    {0.681291281, 0.788964712, 0.990363227},
    {0.722193294, 0.813952739, 0.976574709},
    {0.761464949, 0.834302879, 0.956945269},
    {0.798691636, 0.849786142, 0.931688648},
    {0.833466556, 0.860207984, 0.901068838},
    {0.865395197, 0.86541021, 0.865395561},
    {0.897787179, 0.848937047, 0.820880546},
    {0.924127593, 0.827384882, 0.774508472},
    {0.944468518, 0.800927443, 0.726736146},
    {0.958852946, 0.769767752, 0.678007945},
    {0.96732803, 0.734132809, 0.628751763},
    {0.969954137, 0.694266682, 0.579375448},
    {0.966811177, 0.650421156, 0.530263762},
    {0.958003065, 0.602842431, 0.481775914},
    {0.943660866, 0.551750968, 0.434243684},
    {0.923944917, 0.49730856, 0.387970225},
    {0.89904617, 0.439559467, 0.343229596},
    {0.869186849, 0.378313092, 0.300267182},
    {0.834620542, 0.312874446, 0.259301199},
    {0.795631745, 0.24128379, 0.220525627},
    {0.752534934, 0.157246067, 0.184115123},
    {0.705673158, 0.01555616, 0.150232812}
  };

  for (int i = 0; i < 127; i++) {
    double x = i*(32 - 1e-10)/(127 - 1);
    int j = x; x -= j;
    for (int k = 0; k < 3; k++)
      cmap[i][k] = (1. - x)*basemap[j][k] + x*basemap[j+1][k];
  }
}

void gray (double cmap[127][3])
{
  for (int i = 0; i < 127; i++)
    for (int k = 0; k < 3; k++)
      cmap[i][k] = i/(127 - 1.);
}

void randomap (double cmap[127][3])
{
  srand(0);
  for (int i = 0; i < 127; i++)
    for (int k = 0; k < 3; k++)
      cmap[i][k] = (noise() + 1.)/2.;
}

void blue_white_red (double cmap[127][3])
{
  for (int i = 0; i < (127 + 1)/2; i++) {
    cmap[i][0] = i/((127 - 1)/2.);
    cmap[i][1] = i/((127 - 1)/2.);
    cmap[i][2] = 1.;
  }
  for (int i = 0; i < (127 - 1)/2; i++) {
    cmap[i + (127 + 1)/2][0] = 1.;
    cmap[i + (127 + 1)/2][1] = cmap[(127 - 3)/2 - i][1];
    cmap[i + (127 + 1)/2][2] = cmap[(127 - 3)/2 - i][1];
  }
}





typedef struct {
  unsigned char r, g, b;
} color;

color colormap_color (double cmap[127][3],
        double val, double min, double max)
{
  color c;
  if (val == HUGE) {
    c.r = c.g = c.b = 0;
    return c;
  }
  int i;
  double coef;
  if (max != min)
    val = (val - min)/(max - min);
  else
    val = 0.;
  if (val <= 0.) i = 0, coef = 0.;
  else if (val >= 1.) i = 127 - 2, coef = 1.;
  else {
    i = val*(127 - 1);
    coef = val*(127 - 1) - i;
  }
  if (!(i < 127 - 1)) qassert ("/u/basilisk/src/output.h", 0, "i < NCMAP - 1");
  unsigned char * c1 = (unsigned char *) &c;
  for (int j = 0; j < 3; j++)
    c1[j] = 255*(cmap[i][j]*(1. - coef) + cmap[i + 1][j]*coef);
  return c;
}
#line 340 "/u/basilisk/src/output.h"
static const char * extension (const char * file, const char * ext) {
  int len = strlen(file);
  return len > 4 && !strcmp (file + len - 4, ext) ? file + len - 4 : NULL;
}

static const char * is_animation (const char * file) {
  const char * ext;
  if ((ext = extension (file, ".mp4")) ||
      (ext = extension (file, ".ogv")) ||
      (ext = extension (file, ".gif")))
    return ext;
  return NULL;
}

static struct {
  FILE ** fp;
  char ** names;
  int n;
} open_image_data = {NULL, NULL, 0};

static void open_image_cleanup()
{
  for (int i = 0; i < open_image_data.n; i++) {
    qpclose (open_image_data.fp[i]);
    pfree (open_image_data.names[i],__func__,__FILE__,0);
  }
  pfree (open_image_data.fp,__func__,__FILE__,0);
  pfree (open_image_data.names,__func__,__FILE__,0);
  open_image_data.fp = NULL;
  open_image_data.names = NULL;
  open_image_data.n = 0;
}

static FILE * open_image_lookup (const char * file)
{
  for (int i = 0; i < open_image_data.n; i++)
    if (!strcmp (file, open_image_data.names[i]))
      return open_image_data.fp[i];
  return NULL;
}

static bool which (const char * command)
{
  char * s = getenv ("PATH");
  if (!s)
    return false;
  char path[strlen(s) + 1];
  strcpy (path, s);
  s = strtok (path, ":");
  while (s) {
    char f[strlen(s) + strlen(command) + 2];
    strcpy (f, s);
    strcat (f, "/");
    strcat (f, command);
    FILE * fp = fopen (f, "r");
    if (fp) {
      fclose (fp);
      return true;
    }
    s = strtok (NULL, ":");
  }
  return false;
}

static FILE * ppm_fallback (const char * file, const char * mode)
{
  char filename[strlen(file) + 5];
  strcpy (filename, file);
  strcat (filename, ".ppm");
  FILE * fp = fopen (filename, mode);
  if (!fp) {
    perror (file);



    exit (1);
  }
  return fp;
}

FILE * open_image (const char * file, const char * options)
{
  if (!(pid() == 0)) qassert ("/u/basilisk/src/output.h", 0, "pid() == 0");
  const char * ext;
  if ((ext = is_animation (file))) {
    FILE * fp = open_image_lookup (file);
    if (fp)
      return fp;

    int len = strlen ("ppm2???    ") + strlen (file) +
      (options ? strlen (options) : 0);
    char command[len];
    strcpy (command, "ppm2"); strcat (command, ext + 1);

    static int has_ffmpeg = -1;
    if (has_ffmpeg < 0) {
      if (which (command) && (which ("ffmpeg") || which ("avconv")))
 has_ffmpeg = true;
      else {
 fprintf (ferr,
   "open_image(): cannot find '%s' or 'ffmpeg'/'avconv'\n"
   "  falling back to raw PPM outputs\n", command);
 has_ffmpeg = false;
      }
    }
    if (!has_ffmpeg)
      return ppm_fallback (file, "a");

    static bool added = false;
    if (!added) {
      free_solver_func_add (open_image_cleanup);
      added = true;
    }
    open_image_data.n++;
    open_image_data.names = (char * *) prealloc (open_image_data.names, (open_image_data.n)*sizeof(char *),__func__,__FILE__,0);
    open_image_data.names[open_image_data.n - 1] = pstrdup (file,__func__,__FILE__,0);

    if (options) {
      strcat (command, " ");
      strcat (command, options);
    }
    strcat (command, !strcmp (ext, ".mp4") ? " " : " > ");
    strcat (command, file);
    open_image_data.fp = (FILE * *) prealloc (open_image_data.fp, (open_image_data.n)*sizeof(FILE *),__func__,__FILE__,0);
    return open_image_data.fp[open_image_data.n - 1] = qpopen (command, "w");
  }
  else {
    static int has_convert = -1;
    if (has_convert < 0) {
      if (which ("convert"))
 has_convert = true;
      else {
 fprintf (ferr,
   "open_image(): cannot find 'convert'\n"
   "  falling back to raw PPM outputs\n");
 has_convert = false;
      }
    }
    if (!has_convert)
      return ppm_fallback (file, "w");

    int len = strlen ("convert ppm:-   ") + strlen (file) +
      (options ? strlen (options) : 0);
    char command[len];
    strcpy (command, "convert ppm:- ");
    if (options) {
      strcat (command, options);
      strcat (command, " ");
    }
    strcat (command, file);
    return qpopen (command, "w");
  }
}

void close_image (const char * file, FILE * fp)
{
  if (!(pid() == 0)) qassert ("/u/basilisk/src/output.h", 0, "pid() == 0");
  if (is_animation (file)) {
    if (!open_image_lookup (file))
      fclose (fp);
  }
  else if (which ("convert"))
    qpclose (fp);
  else
    fclose (fp);
}
#line 571 "/u/basilisk/src/output.h"
struct OutputPPM {
  scalar f;
  FILE * fp;
  int n;
  char * file;
  double min, max, spread, z;
  bool linear;
  double box[2][2];
  scalar mask;
  colormap map;
  char * opt;
};

     
void output_ppm (struct OutputPPM p)
{tracing("output_ppm","/u/basilisk/src/output.h",0);

  if (p.n == 0) p.n = N;
  if (p.min == 0 && p.max == 0) {
    stats s = statsf (p.f);
    if (p.spread < 0.)
      p.min = s.min, p.max = s.max;
    else {
      double avg = s.sum/s.volume, spread = (p.spread ? p.spread : 5.)*s.stddev;
      p.min = avg - spread; p.max = avg + spread;
    }
  }
  if (p.box[0][0] == 0. && p.box[0][1] == 0. &&
      p.box[1][0] == 0. && p.box[1][1] == 0.) {
    p.box[0][0] = X0; p.box[0][1] = Y0;
    p.box[1][0] = X0 + L0; p.box[1][1] = Y0 + L0;
  }
  if (!p.map)
    p.map = jet;
  if (p.linear) {
    scalar f = p.f, mask = p.mask;
    if (mask.i)
      boundary_internal ((scalar *)((scalar[]){f, mask,{-1}}), "/u/basilisk/src/output.h", 0);
    else
      boundary_internal ((scalar *)((scalar[]){f,{-1}}), "/u/basilisk/src/output.h", 0);
  }

  double fn = p.n;
  double Delta = (p.box[1][0] - p.box[0][0])/fn;
  int ny = (p.box[1][1] - p.box[0][1])/Delta;
  if (ny % 2) ny++;

  color ** ppm = (color **) matrix_new (ny, p.n, sizeof(color));
  double cmap[127][3];
  p.map (cmap);
  OMP_PARALLEL() {
    OMP(omp for schedule(static))
      for (int j = 0; j < ny; j++) {
 double yp = Delta*j + p.box[0][1] + Delta/2.;
 for (int i = 0; i < p.n; i++) {
   double xp = Delta*i + p.box[0][0] + Delta/2., v;
   if (p.mask.i) {
     if (p.linear) {
       double m = interpolate ((struct _interpolate){p.mask, xp, yp, p.z});
       if (m < 0.)
  v = HUGE;
       else
  v = interpolate ((struct _interpolate){p.f, xp, yp, p.z});
     }
     else {
       Point point = locate ((struct _locate){xp, yp, p.z});int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
       if (point.level < 0 || val(p.mask,0,0,0) < 0.)
  v = HUGE;
       else
  v = val(p.f,0,0,0);
     }
   }
   else if (p.linear)
     v = interpolate ((struct _interpolate){p.f, xp, yp, p.z});
   else {
     Point point = locate ((struct _locate){xp, yp, p.z});int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
     v = point.level >= 0 ? val(p.f,0,0,0) : HUGE;
   }
   ppm[ny - 1 - j][i] = colormap_color (cmap, v, p.min, p.max);
 }
      }
  }

  if (pid() == 0) {
#if _MPI
    MPI_Reduce (MPI_IN_PLACE, ppm[0], 3*ny*p.n, MPI_UNSIGNED_CHAR, MPI_MAX, 0,
  MPI_COMM_WORLD);
#endif
    if (!p.fp) p.fp = fout;
    if (p.file)
      p.fp = open_image (p.file, p.opt);

    fprintf (p.fp, "P6\n%u %u 255\n", p.n, ny);
    fwrite (((void **) ppm)[0], sizeof(color), ny*p.n, p.fp);

    if (p.file)
      close_image (p.file, p.fp);
    else
      fflush (p.fp);
  }
#if _MPI
  else
    MPI_Reduce (ppm[0], NULL, 3*ny*p.n, MPI_UNSIGNED_CHAR, MPI_MAX, 0,
  MPI_COMM_WORLD);
#endif

  matrix_free (ppm);
end_tracing("output_ppm","/u/basilisk/src/output.h",0);}
#line 710 "/u/basilisk/src/output.h"
struct OutputGRD {
  scalar f;
  FILE * fp;
  double Delta;
  bool linear;
  double box[2][2];
  scalar mask;
};

     
void output_grd (struct OutputGRD p)
{tracing("output_grd","/u/basilisk/src/output.h",0);

  if (!p.fp) p.fp = fout;
  if (p.box[0][0] == 0. && p.box[0][1] == 0. &&
      p.box[1][0] == 0. && p.box[1][1] == 0.) {
    p.box[0][0] = X0; p.box[0][1] = Y0;
    p.box[1][0] = X0 + L0; p.box[1][1] = Y0 + L0;
    if (p.Delta == 0) p.Delta = L0/N;
  }
  if (p.linear) {
    scalar f = p.f, mask = p.mask;
    if (mask.i)
      boundary_internal ((scalar *)((scalar[]){f, mask,{-1}}), "/u/basilisk/src/output.h", 0);
    else
      boundary_internal ((scalar *)((scalar[]){f,{-1}}), "/u/basilisk/src/output.h", 0);
  }

  double Delta = p.Delta;
  int nx = (p.box[1][0] - p.box[0][0])/Delta;
  int ny = (p.box[1][1] - p.box[0][1])/Delta;


  fprintf (p.fp, "ncols          %d\n", nx);
  fprintf (p.fp, "nrows          %d\n", ny);
  fprintf (p.fp, "xllcorner      %g\n", p.box[0][0]);
  fprintf (p.fp, "yllcorner      %g\n", p.box[0][1]);
  fprintf (p.fp, "cellsize       %g\n", Delta);
  fprintf (p.fp, "nodata_value   -9999\n");


  for (int j = ny-1; j >= 0; j--) {
    double yp = Delta*j + p.box[0][1] + Delta/2.;
    for (int i = 0; i < nx; i++) {
      double xp = Delta*i + p.box[0][0] + Delta/2., v;
      if (p.mask.i) {
 if (p.linear) {
   double m = interpolate ((struct _interpolate){p.mask, xp, yp});
   if (m < 0.)
     v = HUGE;
   else
     v = interpolate ((struct _interpolate){p.f, xp, yp});
 }
 else {
   Point point = locate ((struct _locate){xp, yp});int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
   if (point.level < 0 || val(p.mask,0,0,0) < 0.)
     v = HUGE;
   else
     v = val(p.f,0,0,0);
 }
      }
      else if (p.linear)
 v = interpolate ((struct _interpolate){p.f, xp, yp});
      else {
 Point point = locate ((struct _locate){xp, yp});int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
 v = point.level >= 0 ? val(p.f,0,0,0) : HUGE;
      }
      if (v == HUGE)
 fprintf (p.fp, "-9999 ");
      else
 fprintf (p.fp, "%f ", v);
    }
    fprintf (p.fp, "\n");
  }

  fflush (p.fp);
end_tracing("output_grd","/u/basilisk/src/output.h",0);}
#line 813 "/u/basilisk/src/output.h"
struct OutputGfs {
  FILE * fp;
  scalar * list;
  double t;
  char * file;
  bool translate;
};

static char * replace (const char * input, int target, int with,
         bool translate)
{
  if (translate) {
    if (!strcmp (input, "u.x"))
      return pstrdup ("U",__func__,__FILE__,0);
    if (!strcmp (input, "u.y"))
      return pstrdup ("V",__func__,__FILE__,0);
    if (!strcmp (input, "u.z"))
      return pstrdup ("W",__func__,__FILE__,0);
  }
  char * name = pstrdup (input,__func__,__FILE__,0), * i = name;
  while (*i != '\0') {
    if (*i == target)
      *i = with;
    i++;
  }
  return name;
}

     
void output_gfs (struct OutputGfs p)
{tracing("output_gfs","/u/basilisk/src/output.h",0);
  char * fname = p.file;

#if _MPI



  FILE * fp = p.fp;
  if (p.file == NULL) {
    long pid = getpid();
    MPI_Bcast (&pid, 1, MPI_LONG, 0, MPI_COMM_WORLD);
    fname = ((char *) pmalloc ((80)*sizeof(char),__func__,__FILE__,0));
    snprintf (fname, 80, ".output-%ld", pid);
    p.fp = NULL;
  }
#endif

  bool opened = false;
  if (p.fp == NULL) {
    if (fname == NULL)
      p.fp = fout;
    else if (!(p.fp = fopen (fname, "w"))) {
      perror (fname);
      exit (1);
    }
    else
      opened = true;
  }

  scalar * list = p.list ? p.list : list_copy (all);

  restriction (list);
  fprintf (p.fp,
    "1 0 GfsSimulation GfsBox GfsGEdge { binary = 1"
    " x = %g y = %g ",
    0.5 + X0/L0, 0.5 + Y0/L0);




  if (list != NULL && list[0].i != -1) {
    scalar s = list[0];
    char * name = replace (_attribute[s.i].name, '.', '_', p.translate);
    fprintf (p.fp, "variables = %s", name);
    pfree (name,__func__,__FILE__,0);
    for (int i = 1; i < list_len(list); i++) {
      scalar s = list[i];
      if (_attribute[s.i].name) {
 char * name = replace (_attribute[s.i].name, '.', '_', p.translate);
 fprintf (p.fp, ",%s", name);
 pfree (name,__func__,__FILE__,0);
      }
    }
    fprintf (p.fp, " ");
  }
  fprintf (p.fp, "} {\n");
  fprintf (p.fp, "  Time { t = %g }\n", t);
  if (L0 != 1.)
    fprintf (p.fp, "  PhysicalParams { L = %g }\n", L0);
  fprintf (p.fp, "  VariableTracerVOF f\n");
  fprintf (p.fp, "}\nGfsBox { x = 0 y = 0 z = 0 } {\n");

#if _MPI
  long header;
  if ((header = ftell (p.fp)) < 0) {
    perror ("output_gfs(): error in header");
    exit (1);
  }
  int cell_size = sizeof(unsigned) + sizeof(double);
  {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    if (_attribute[s.i].name)
      cell_size += sizeof(double);}}
  scalar index = new_scalar("index");
  size_t total_size = header + (z_indexing (index, false) + 1)*cell_size;
#endif



  {foreach_cell() {
#if _MPI
    if (is_local(cell))
#endif
    {
#if _MPI
      if (fseek (p.fp, header + val(index,0,0,0)*cell_size, SEEK_SET) < 0) {
 perror ("output_gfs(): error while seeking");
 exit (1);
      }
#endif
      unsigned flags =
 level == 0 ? 0 :



      child.x == -1 && child.y == -1 ? 0 :
 child.x == -1 && child.y == 1 ? 1 :
 child.x == 1 && child.y == -1 ? 2 :
 3;
#line 951 "/u/basilisk/src/output.h"
      if (is_leaf(cell))
 flags |= (1 << 4);
      fwrite (&flags, sizeof (unsigned), 1, p.fp);
      double a = -1;
      fwrite (&a, sizeof (double), 1, p.fp);
      {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
 if (_attribute[s.i].name) {
   if (_attribute[s.i].v.x.i >= 0) {




     if (_attribute[s.i].v.x.i == s.i) {
       s = _attribute[s.i].v.y;
       a = is_local(cell) && val(s,0,0,0) != HUGE ? val(s,0,0,0) : (double) DBL_MAX;
     }
     else if (_attribute[s.i].v.y.i == s.i) {
       s = _attribute[s.i].v.x;
       a = is_local(cell) && val(s,0,0,0) != HUGE ? - val(s,0,0,0) : (double) DBL_MAX;
     }





   }
   else
     a = is_local(cell) && val(s,0,0,0) != HUGE ? val(s,0,0,0) : (double) DBL_MAX;
   fwrite (&a, sizeof (double), 1, p.fp);
 }}}
    }
    if (is_leaf(cell))
      continue;
  }end_foreach_cell();}

#if _MPI
  delete (((scalar[]){index,{-1}}));
  if (!pid() && fseek (p.fp, total_size, SEEK_SET) < 0) {
    perror ("output_gfs(): error while finishing");
    exit (1);
  }
  if (!pid())
#endif
    fputs ("}\n", p.fp);
  fflush (p.fp);

  if (!p.list)
    pfree (list,__func__,__FILE__,0);
  if (opened)
    fclose (p.fp);

#if _MPI
  if (p.file == NULL) {
    MPI_Barrier (MPI_COMM_WORLD);
    if (pid() == 0) {
      if (fp == NULL)
 fp = fout;
      p.fp = fopen (fname, "r");
      size_t l;
      unsigned char buffer[8192];
      while ((l = fread (buffer, 1, 8192, p.fp)) > 0)
 fwrite (buffer, 1, l, fp);
      fflush (fp);
      remove (fname);
    }
    pfree (fname,__func__,__FILE__,0);
  }
#endif
end_tracing("output_gfs","/u/basilisk/src/output.h",0);}
#line 1043 "/u/basilisk/src/output.h"
struct Dump {
  char * file;
  scalar * list;
  FILE * fp;
  bool unbuffered;
};

struct DumpHeader {
  double t;
  long len;
  int i, depth, npe, version;
  coord n;
};

static const int dump_version =

  170901;

static scalar * dump_list (scalar * lista)
{
  scalar * list = is_constant(cm) ? NULL : list_concat (((scalar[]){cm,{-1}}), NULL);
  {scalar*_i=(scalar*)( lista);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    if (!_attribute[s.i].face && !_attribute[s.i].nodump && s.i != cm.i)
      list = list_add (list, s);}}
  return list;
}

static void dump_header (FILE * fp, struct DumpHeader * header, scalar * list)
{
  if (fwrite (header, sizeof(struct DumpHeader), 1, fp) < 1) {
    perror ("dump(): error while writing header");
    exit (1);
  }
  {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {
    unsigned len = strlen(_attribute[s.i].name);
    if (fwrite (&len, sizeof(unsigned), 1, fp) < 1) {
      perror ("dump(): error while writing len");
      exit (1);
    }
    if (fwrite (_attribute[s.i].name, sizeof(char), len, fp) < len) {
      perror ("dump(): error while writing s.name");
      exit (1);
    }
  }}}
  double o[4] = {X0,Y0,Z0,L0};
  if (fwrite (o, sizeof(double), 4, fp) < 4) {
    perror ("dump(): error while writing coordinates");
    exit (1);
  }
}

#if !_MPI
     
void dump (struct Dump p)
{tracing("dump","/u/basilisk/src/output.h",0);
  FILE * fp = p.fp;
  char def[] = "dump", * file = p.file ? p.file : p.fp ? NULL : def;

  char * name = NULL;
  if (file) {
    name = (char *) pmalloc (strlen(file) + 2,__func__,__FILE__,0);
    strcpy (name, file);
    if (!p.unbuffered)
      strcat (name, "~");
    if ((fp = fopen (name, "w")) == NULL) {
      perror (name);
      exit (1);
    }
  }
  if (!(fp)) qassert ("/u/basilisk/src/output.h", 0, "fp");

  scalar * dlist = dump_list (p.list ? p.list : all);
  scalar  size=new_scalar("size");
  scalar * list = list_concat (((scalar[]){size,{-1}}), dlist); pfree (dlist,__func__,__FILE__,0);
  struct DumpHeader header = { t, list_len(list), iter, depth(), npe(),
          dump_version };
  dump_header (fp, &header, list);

  subtree_size (size, false);

  {foreach_cell() {
    unsigned flags = is_leaf(cell) ? leaf : 0;
    if (fwrite (&flags, sizeof(unsigned), 1, fp) < 1) {
      perror ("dump(): error while writing flags");
      exit (1);
    }
    {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      if (fwrite (&val(s,0,0,0), sizeof(double), 1, fp) < 1) {
 perror ("dump(): error while writing scalars");
 exit (1);
      }}}
    if (is_leaf(cell))
      continue;
  }end_foreach_cell();}

  pfree (list,__func__,__FILE__,0);
  if (file) {
    fclose (fp);
    if (!p.unbuffered)
      rename (name, file);
    pfree (name,__func__,__FILE__,0);
  }delete((scalar*)((scalar[]){size,{-1}}));
end_tracing("dump","/u/basilisk/src/output.h",0);}
#else
     
void dump (struct Dump p)
{tracing("dump","/u/basilisk/src/output.h",0);
  FILE * fp = p.fp;
  char def[] = "dump", * file = p.file ? p.file : p.fp ? NULL : def;

  if (fp != NULL || file == NULL) {
    fprintf (ferr, "dump(): must specify a file name when using MPI\n");
    exit(1);
  }

  char name[strlen(file) + 2];
  strcpy (name, file);
  if (!p.unbuffered)
    strcat (name, "~");
  FILE * fh = fopen (name, "w");
  if (fh == NULL) {
    perror (name);
    exit (1);
  }

  scalar * dlist = dump_list (p.list ? p.list : all);
  scalar  size=new_scalar("size");
  scalar * list = list_concat (((scalar[]){size,{-1}}), dlist); pfree (dlist,__func__,__FILE__,0);
  struct DumpHeader header = { t, list_len(list), iter, depth(), npe(),
          dump_version };







  if (pid() == 0)
    dump_header (fh, &header, list);

  scalar index = {-1};

  index = new_scalar("index");
  z_indexing (index, false);
  int cell_size = sizeof(unsigned) + header.len*sizeof(double);
  int sizeofheader = sizeof(header) + 4*sizeof(double);
  {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    sizeofheader += sizeof(unsigned) + sizeof(char)*strlen(_attribute[s.i].name);}}
  long pos = pid() ? 0 : sizeofheader;

  subtree_size (size, false);

  {foreach_cell() {

    if (is_local(cell)) {
      long offset = sizeofheader + val(index,0,0,0)*cell_size;
      if (pos != offset) {
 fseek (fh, offset, SEEK_SET);
 pos = offset;
      }
      unsigned flags = is_leaf(cell) ? leaf : 0;
      fwrite (&flags, 1, sizeof(unsigned), fh);
      {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
 fwrite (&val(s,0,0,0), 1, sizeof(double), fh);}}
      pos += cell_size;
    }
    if (is_leaf(cell))
      continue;
  }end_foreach_cell();}

  delete (((scalar[]){index,{-1}}));

  pfree (list,__func__,__FILE__,0);
  fclose (fh);
  if (!p.unbuffered && pid() == 0)
    rename (name, file);delete((scalar*)((scalar[]){size,{-1}}));
end_tracing("dump","/u/basilisk/src/output.h",0);}
#endif

     
bool restore (struct Dump p)
{tracing("restore","/u/basilisk/src/output.h",0);
  FILE * fp = p.fp;
  char * file = p.file;
  if (file && (fp = fopen (file, "r")) == NULL)
    {end_tracing("restore","/u/basilisk/src/output.h",0);return false;}
  if (!(fp)) qassert ("/u/basilisk/src/output.h", 0, "fp");

  struct DumpHeader header;
  if (fread (&header, sizeof(header), 1, fp) < 1) {
    fprintf (ferr, "restore(): error: expecting header\n");
    exit (1);
  }


  init_grid (1);
  {foreach_cell() {
    cell.pid = pid();
    cell.flags |= active;
  }end_foreach_cell();}
  ((Tree *)grid)->dirty = true;
#line 1264 "/u/basilisk/src/output.h"
  bool restore_all = (p.list == all);
  scalar * list = dump_list (p.list ? p.list : all);
  if (header.version == 161020) {
    if (header.len - 1 != list_len (list)) {
      fprintf (ferr,
        "restore(): error: the list lengths don't match: "
        "%ld (file) != %d (code)\n",
        header.len - 1, list_len (list));
      exit (1);
    }
  }
  else {
    if (header.version != dump_version) {
      fprintf (ferr,
        "restore(): error: file version mismatch: "
        "%d (file) != %d (code)\n",
        header.version, dump_version);
      exit (1);
    }

    scalar * input = NULL;
    for (int i = 0; i < header.len; i++) {
      unsigned len;
      if (fread (&len, sizeof(unsigned), 1, fp) < 1) {
 fprintf (ferr, "restore(): error: expecting len\n");
 exit (1);
      }
      char name[len + 1];
      if (fread (name, sizeof(char), len, fp) < 1) {
 fprintf (ferr, "restore(): error: expecting s.name\n");
 exit (1);
      }
      name[len] = '\0';

      if (i > 0) {
 bool found = false;
 {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
   if (!strcmp (_attribute[s.i].name, name)) {
     input = list_append (input, s);
     found = true; break;
   }}}
 if (!found) {
   if (restore_all) {
     scalar s = new_scalar("s");
     pfree (_attribute[s.i].name,__func__,__FILE__,0);
     _attribute[s.i].name = pstrdup (name,__func__,__FILE__,0);
     input = list_append (input, s);
   }
   else
     input = list_append (input, (scalar){INT_MAX});
 }
      }
    }
    pfree (list,__func__,__FILE__,0);
    list = input;

    double o[4];
    if (fread (o, sizeof(double), 4, fp) < 4) {
      fprintf (ferr, "restore(): error: expecting coordinates\n");
      exit (1);
    }
    origin ((struct _origin){o[0], o[1], o[2]});
    size (o[3]);
  }
#line 1339 "/u/basilisk/src/output.h"
  scalar * listm = is_constant(cm) ? NULL : (scalar *)((vector[]){fm,{{-1},{-1}}});



  {foreach_cell() {
    unsigned flags;
    if (fread (&flags, sizeof(unsigned), 1, fp) != 1) {
      fprintf (ferr, "restore(): error: expecting 'flags'\n");
      exit (1);
    }

    fseek (fp, sizeof(double), SEEK_CUR);
    {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {
      double val;
      if (fread (&val, sizeof(double), 1, fp) != 1) {
 fprintf (ferr, "restore(): error: expecting a scalar\n");
 exit (1);
      }
      if (s.i != INT_MAX)
 val(s,0,0,0) = val;
    }}}
    if (!(flags & leaf) && is_leaf(cell))
      refine_cell (point, listm, 0, NULL);
    if (is_leaf(cell))
      continue;
  }end_foreach_cell();}
  {scalar*_i=(scalar*)( all);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    _attribute[s.i].dirty = true;}}


  scalar * other = NULL;
  {scalar*_i=(scalar*)( all);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    if (!list_lookup (list, s) && !list_lookup (listm, s))
      other = list_append (other, s);}}
  reset (other, 0.);
  pfree (other,__func__,__FILE__,0);

  pfree (list,__func__,__FILE__,0);
  if (file)
    fclose (fp);


  while (iter < header.i && events (false))
    iter = inext;
  events (false);
  while (t < header.t && events (false))
    t = tnext;
  t = header.t;
  events (false);

  {end_tracing("restore","/u/basilisk/src/output.h",0);return true;}
end_tracing("restore","/u/basilisk/src/output.h",0);}
#line 431 "/u/basilisk/src/utils.h"
#line 12 "/u/basilisk/src/run.h"

     
void run (void)
{tracing("run","/u/basilisk/src/run.h",0);
  iter = 0, t = 0., dt = 1.;
  init_grid (N);

  perf.nc = perf.tnc = 0;
  perf.gt = timer_start();
  while (events (true)) {





    update_perf();
    iter = inext, t = tnext;
  }




  timer_print (perf.gt, iter, perf.tnc);

  free_grid();
end_tracing("run","/u/basilisk/src/run.h",0);}




static int defaults_0_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i = 0);*ip=i;*tp=t;return ret;}      static int defaults_0(const int i,const double t,Event *_ev){tracing("defaults_0","/u/basilisk/src/run.h",0); {
  display ((struct _display){"box();"});
}{end_tracing("defaults_0","/u/basilisk/src/run.h",0);return 0;}end_tracing("defaults_0","/u/basilisk/src/run.h",0);}





static int cleanup_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(t = 1234567890);*ip=i;*tp=t;return ret;}      static int cleanup(const int i,const double t,Event *_ev){tracing("cleanup","/u/basilisk/src/run.h",0); {
  display ((struct _display){"", true});
}{end_tracing("cleanup","/u/basilisk/src/run.h",0);return 0;}end_tracing("cleanup","/u/basilisk/src/run.h",0);}
#line 28 "/u/basilisk/src/navier-stokes/centered.h"
#line 1 "./timestep.h"
#line 1 "/u/basilisk/src/timestep.h"

double timestep (const vector u, double dtmax)
{
  static double previous = 0.;
  dtmax /= CFL;
  foreach_face_stencil(){_stencil_is_face_x(){
    {_stencil_val(u.x,0,0,0); {   
      _stencil_val(u.x,0,0,0);

_stencil_val(fm.x,0,0,0); 
_stencil_val(fm.x,0,0,0);    

          
       



         
    
#line 16
}   }}end__stencil_is_face_x()
#line 6
_stencil_is_face_y(){
    {_stencil_val(u.y,0,0,0); {   
      _stencil_val(u.y,0,0,0);

_stencil_val(fm.y,0,0,0); 
_stencil_val(fm.y,0,0,0);    

          
       



         
    
#line 16
}   }}end__stencil_is_face_y()}end_foreach_face_stencil();
  
#line 6
if(!is_constant(fm.x)){
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel reduction(min:dtmax)){
#line 6
foreach_face_generic(){is_face_x(){
    if (val(u.x,0,0,0) != 0.) {
      double dt = Delta/fabs(val(u.x,0,0,0));

      if (!(val(fm.x,0,0,0))) qassert ("/u/basilisk/src/timestep.h", 0, "fm.x[]");
      dt *= val(fm.x,0,0,0);



      if (dt < dtmax) dtmax = dt;
    }}end_is_face_x()
#line 6
is_face_y(){
    if (val(u.y,0,0,0) != 0.) {
      double dt = Delta/fabs(val(u.y,0,0,0));

      if (!(val(fm.y,0,0,0))) qassert ("/u/basilisk/src/timestep.h", 0, "fm.x[]");
      dt *= val(fm.y,0,0,0);



      if (dt < dtmax) dtmax = dt;
    }}end_is_face_y()}end_foreach_face_generic();mpi_all_reduce_array(&dtmax,double,MPI_MIN,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 16
}else {struct{double x,y;}_const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);
  
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel reduction(min:dtmax)){
#line 6
foreach_face_generic(){is_face_x(){
    if (val(u.x,0,0,0) != 0.) {
      double dt = Delta/fabs(val(u.x,0,0,0));

      if (!(_const_fm.x)) qassert ("/u/basilisk/src/timestep.h", 0, "fm.x[]");
      dt *= _const_fm.x;



      if (dt < dtmax) dtmax = dt;
    }}end_is_face_x()
#line 6
is_face_y(){
    if (val(u.y,0,0,0) != 0.) {
      double dt = Delta/fabs(val(u.y,0,0,0));

      if (!(_const_fm.y)) qassert ("/u/basilisk/src/timestep.h", 0, "fm.x[]");
      dt *= _const_fm.y;



      if (dt < dtmax) dtmax = dt;
    }}end_is_face_y()}end_foreach_face_generic();mpi_all_reduce_array(&dtmax,double,MPI_MIN,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 16
}
  dtmax *= CFL;
  if (dtmax > previous)
    dtmax = (previous + 0.1*dtmax)/1.1;
  previous = dtmax;
  return dtmax;
}
#line 29 "/u/basilisk/src/navier-stokes/centered.h"
#line 1 "./bcg.h"
#line 1 "/u/basilisk/src/bcg.h"
#line 11 "/u/basilisk/src/bcg.h"
void tracer_fluxes (scalar f,
      vector uf,
      vector flux,
      double dt,
              scalar src)
{





  vector  g=new_vector("g");
  gradients (((scalar[]){f,{-1}}), ((vector[]){g,{{-1},{-1}}}));




  foreach_face_stencil(){_stencil_is_face_x(){ {        







    _stencil_val(fm.x,0,0,0);_stencil_val(uf.x,0,0,0);              
    
    _stencil_val(g.x,o_stencil,0,0); _stencil_val(src,-1,0,0);_stencil_val(src,0,0,0);_stencil_val(f, o_stencil,0,0);





_stencil_val(fm.y,o_stencil,0,0); _stencil_val(fm.y,o_stencil,1,0); {     
       _stencil_val(fm.y,o_stencil,1,0);_stencil_val(fm.y,o_stencil,0,0); _stencil_val(uf.y,o_stencil,1,0);_stencil_val(uf.y,o_stencil,0,0);         
       _stencil_val(f,o_stencil,-1,0);_stencil_val(f, o_stencil,0,0);_stencil_val(f, o_stencil,0,0); _stencil_val(f,o_stencil,1,0);
        
    }





      
#line 58 "/u/basilisk/src/bcg.h"
    _stencil_val_a(flux.x,0,0,0);_stencil_val(uf.x,0,0,0);  
  }}end__stencil_is_face_x()
#line 28
_stencil_is_face_y(){ {        







    _stencil_val(fm.y,0,0,0);_stencil_val(uf.y,0,0,0);              
    
    _stencil_val(g.y,0,o_stencil,0); _stencil_val(src,0,-1,0);_stencil_val(src,0,0,0);_stencil_val(f,0, o_stencil,0);





_stencil_val(fm.x,0,o_stencil,0); _stencil_val(fm.x,1,o_stencil,0); {     
       _stencil_val(fm.x,1,o_stencil,0);_stencil_val(fm.x,0,o_stencil,0); _stencil_val(uf.x,1,o_stencil,0);_stencil_val(uf.x,0,o_stencil,0);         
       _stencil_val(f,-1,o_stencil,0);_stencil_val(f,0, o_stencil,0);_stencil_val(f,0, o_stencil,0); _stencil_val(f,1,o_stencil,0);
        
    }





      
#line 58 "/u/basilisk/src/bcg.h"
    _stencil_val_a(flux.y,0,0,0);_stencil_val(uf.y,0,0,0);  
  }}end__stencil_is_face_y()}end_foreach_face_stencil();




  
#line 28
if(!is_constant(fm.x) && !is_constant(src)){{foreach_face_generic(){is_face_x(){ {







    double un = dt*val(uf.x,0,0,0)/(val(fm.x,0,0,0)*Delta + 1e-30), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = val(f,i,0,0) + (val(src,0,0,0) + val(src,-1,0,0))*dt/4. + s*(1. - s*un)*val(g.x,i,0,0)*Delta/2.;





    if (val(fm.y,i,0,0) && val(fm.y,i,1,0)) {
      double vn = (val(uf.y,i,0,0) + val(uf.y,i,1,0))/(val(fm.y,i,0,0) + val(fm.y,i,1,0));
      double fyy = vn < 0. ? val(f,i,1,0) - val(f,i,0,0) : val(f,i,0,0) - val(f,i,-1,0);
      f2 -= dt*vn*fyy/(2.*Delta);
    }
#line 58 "/u/basilisk/src/bcg.h"
    val(flux.x,0,0,0) = f2*val(uf.x,0,0,0);
  }}end_is_face_x()
#line 28
is_face_y(){ {







    double un = dt*val(uf.y,0,0,0)/(val(fm.y,0,0,0)*Delta + 1e-30), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = val(f,0,i,0) + (val(src,0,0,0) + val(src,0,-1,0))*dt/4. + s*(1. - s*un)*val(g.y,0,i,0)*Delta/2.;





    if (val(fm.x,0,i,0) && val(fm.x,1,i,0)) {
      double vn = (val(uf.x,0,i,0) + val(uf.x,1,i,0))/(val(fm.x,0,i,0) + val(fm.x,1,i,0));
      double fyy = vn < 0. ? val(f,1,i,0) - val(f,0,i,0) : val(f,0,i,0) - val(f,-1,i,0);
      f2 -= dt*vn*fyy/(2.*Delta);
    }
#line 58 "/u/basilisk/src/bcg.h"
    val(flux.y,0,0,0) = f2*val(uf.y,0,0,0);
  }}end_is_face_y()}end_foreach_face_generic();}}else if(is_constant(fm.x) && !is_constant(src)){struct{double x,y;}_const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);




  {
#line 28
foreach_face_generic(){is_face_x(){ {







    double un = dt*val(uf.x,0,0,0)/(_const_fm.x*Delta + 1e-30), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = val(f,i,0,0) + (val(src,0,0,0) + val(src,-1,0,0))*dt/4. + s*(1. - s*un)*val(g.x,i,0,0)*Delta/2.;





    if (_const_fm.y && _const_fm.y) {
      double vn = (val(uf.y,i,0,0) + val(uf.y,i,1,0))/(_const_fm.y + _const_fm.y);
      double fyy = vn < 0. ? val(f,i,1,0) - val(f,i,0,0) : val(f,i,0,0) - val(f,i,-1,0);
      f2 -= dt*vn*fyy/(2.*Delta);
    }
#line 58 "/u/basilisk/src/bcg.h"
    val(flux.x,0,0,0) = f2*val(uf.x,0,0,0);
  }}end_is_face_x()
#line 28
is_face_y(){ {







    double un = dt*val(uf.y,0,0,0)/(_const_fm.y*Delta + 1e-30), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = val(f,0,i,0) + (val(src,0,0,0) + val(src,0,-1,0))*dt/4. + s*(1. - s*un)*val(g.y,0,i,0)*Delta/2.;





    if (_const_fm.x && _const_fm.x) {
      double vn = (val(uf.x,0,i,0) + val(uf.x,1,i,0))/(_const_fm.x + _const_fm.x);
      double fyy = vn < 0. ? val(f,1,i,0) - val(f,0,i,0) : val(f,0,i,0) - val(f,-1,i,0);
      f2 -= dt*vn*fyy/(2.*Delta);
    }
#line 58 "/u/basilisk/src/bcg.h"
    val(flux.y,0,0,0) = f2*val(uf.y,0,0,0);
  }}end_is_face_y()}end_foreach_face_generic();}}else if(!is_constant(fm.x) && is_constant(src)){double _const_src=_constant[src.i-_NVARMAX];NOT_UNUSED(_const_src);




  {
#line 28
foreach_face_generic(){is_face_x(){ {







    double un = dt*val(uf.x,0,0,0)/(val(fm.x,0,0,0)*Delta + 1e-30), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = val(f,i,0,0) + (_const_src + _const_src)*dt/4. + s*(1. - s*un)*val(g.x,i,0,0)*Delta/2.;





    if (val(fm.y,i,0,0) && val(fm.y,i,1,0)) {
      double vn = (val(uf.y,i,0,0) + val(uf.y,i,1,0))/(val(fm.y,i,0,0) + val(fm.y,i,1,0));
      double fyy = vn < 0. ? val(f,i,1,0) - val(f,i,0,0) : val(f,i,0,0) - val(f,i,-1,0);
      f2 -= dt*vn*fyy/(2.*Delta);
    }
#line 58 "/u/basilisk/src/bcg.h"
    val(flux.x,0,0,0) = f2*val(uf.x,0,0,0);
  }}end_is_face_x()
#line 28
is_face_y(){ {







    double un = dt*val(uf.y,0,0,0)/(val(fm.y,0,0,0)*Delta + 1e-30), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = val(f,0,i,0) + (_const_src + _const_src)*dt/4. + s*(1. - s*un)*val(g.y,0,i,0)*Delta/2.;





    if (val(fm.x,0,i,0) && val(fm.x,1,i,0)) {
      double vn = (val(uf.x,0,i,0) + val(uf.x,1,i,0))/(val(fm.x,0,i,0) + val(fm.x,1,i,0));
      double fyy = vn < 0. ? val(f,1,i,0) - val(f,0,i,0) : val(f,0,i,0) - val(f,-1,i,0);
      f2 -= dt*vn*fyy/(2.*Delta);
    }
#line 58 "/u/basilisk/src/bcg.h"
    val(flux.y,0,0,0) = f2*val(uf.y,0,0,0);
  }}end_is_face_y()}end_foreach_face_generic();}}else {struct{double x,y;}_const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_src=_constant[src.i-_NVARMAX];NOT_UNUSED(_const_src);




  {
#line 28
foreach_face_generic(){is_face_x(){ {







    double un = dt*val(uf.x,0,0,0)/(_const_fm.x*Delta + 1e-30), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = val(f,i,0,0) + (_const_src + _const_src)*dt/4. + s*(1. - s*un)*val(g.x,i,0,0)*Delta/2.;





    if (_const_fm.y && _const_fm.y) {
      double vn = (val(uf.y,i,0,0) + val(uf.y,i,1,0))/(_const_fm.y + _const_fm.y);
      double fyy = vn < 0. ? val(f,i,1,0) - val(f,i,0,0) : val(f,i,0,0) - val(f,i,-1,0);
      f2 -= dt*vn*fyy/(2.*Delta);
    }
#line 58 "/u/basilisk/src/bcg.h"
    val(flux.x,0,0,0) = f2*val(uf.x,0,0,0);
  }}end_is_face_x()
#line 28
is_face_y(){ {







    double un = dt*val(uf.y,0,0,0)/(_const_fm.y*Delta + 1e-30), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = val(f,0,i,0) + (_const_src + _const_src)*dt/4. + s*(1. - s*un)*val(g.y,0,i,0)*Delta/2.;





    if (_const_fm.x && _const_fm.x) {
      double vn = (val(uf.x,0,i,0) + val(uf.x,1,i,0))/(_const_fm.x + _const_fm.x);
      double fyy = vn < 0. ? val(f,1,i,0) - val(f,0,i,0) : val(f,0,i,0) - val(f,-1,i,0);
      f2 -= dt*vn*fyy/(2.*Delta);
    }
#line 58 "/u/basilisk/src/bcg.h"
    val(flux.y,0,0,0) = f2*val(uf.y,0,0,0);
  }}end_is_face_y()}end_foreach_face_generic();}}delete((scalar*)((vector[]){g,{{-1},{-1}}}));
}






struct Advection {
  scalar * tracers;
  vector u;
  double dt;
  scalar * src;
};

void advection (struct Advection p)
{




  scalar * lsrc = p.src;
  if (!lsrc)
    {scalar*_i=(scalar*)( p.tracers);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      lsrc = list_append (lsrc, zeroc);}}
  if (!(list_len(p.tracers) == list_len(lsrc))) qassert ("/u/basilisk/src/bcg.h", 0, "list_len(p.tracers) == list_len(lsrc)");

  scalar f, src;
  {scalar*_i0=lsrc;scalar*_i1= p.tracers;if(_i0)for(src=*_i0,f=*_i1;_i0->i>= 0;src=*++_i0,f=*++_i1){ {
    vector  flux=new_face_vector("flux");
    tracer_fluxes (f, p.u, flux, p.dt, src);





    update_tracer (f, p.u, flux, p.dt);delete((scalar*)((vector[]){flux,{{-1},{-1}}}));

  }}}

  if (!p.src)
    pfree (lsrc,__func__,__FILE__,0);
}
#line 30 "/u/basilisk/src/navier-stokes/centered.h"

#line 1 "./viscosity-embed.h"
#line 1 "/u/basilisk/src/viscosity-embed.h"
#line 1 "./poisson.h"
#line 1 "/u/basilisk/src/poisson.h"
#line 32 "/u/basilisk/src/poisson.h"
void mg_cycle (scalar * a, scalar * res, scalar * da,
        void (* relax) (scalar * da, scalar * res,
          int depth, void * data),
        void * data,
        int nrelax, int minlevel, int maxlevel)
{




  restriction (res);





  minlevel = min (minlevel, maxlevel);
  for (int l = minlevel; l <= maxlevel; l++) {




    if (l == minlevel)
      {foreach_level_or_leaf (l)
 {scalar*_i=(scalar*)( da);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
  
     val(s,0,0,0) = 0.;}}end_foreach_level_or_leaf();}





    else
      {foreach_level (l)
 {scalar*_i=(scalar*)( da);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
  
     val(s,0,0,0) = bilinear_embed(point, s);}}end_foreach_level();}





    boundary_level (da, l);
    for (int i = 0; i < nrelax; i++) {
      relax (da, res, l, data);
      boundary_level (da, l);
    }
  }




  foreach_stencil() {
    scalar s, ds;
    {scalar*_i0= da;scalar*_i1= a;if(_i0)for(ds=*_i0,s=*_i1;_i0->i>= 0;ds=*++_i0,s=*++_i1){
     
 {_stencil_val_r(s,0,0,0); _stencil_val(ds,0,0,0); }}}
  }end_foreach_stencil();




  {
#line 84
foreach() {
    scalar s, ds;
    {scalar*_i0= da;scalar*_i1= a;if(_i0)for(ds=*_i0,s=*_i1;_i0->i>= 0;ds=*++_i0,s=*++_i1){
     
 val(s,0,0,0) += val(ds,0,0,0);}}
  }end_foreach();}
}
#line 102 "/u/basilisk/src/poisson.h"
int NITERMAX = 100, NITERMIN = 1;
double TOLERANCE = 1e-3;




typedef struct {
  int i;
  double resb, resa;
  double sum;
  int nrelax;
  int minlevel;
} mgstats;
#line 125 "/u/basilisk/src/poisson.h"
struct MGSolve {
  scalar * a, * b;
  double (* residual) (scalar * a, scalar * b, scalar * res,
         void * data);
  void (* relax) (scalar * da, scalar * res, int depth,
    void * data);
  void * data;

  int nrelax;
  scalar * res;
  int minlevel;
  double tolerance;
};

mgstats mg_solve (struct MGSolve p)
{





  scalar * da = list_clone (p.a), * res = p.res;
  if (!res)
    res = list_clone (p.b);






  for (int b = 0; b < nboundary; b++)
    {scalar*_i=(scalar*)( da);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      _attribute[s.i].boundary[b] = _attribute[s.i].boundary_homogeneous[b];}}




  mgstats s = {0};
  double sum = 0.;
  foreach_stencil ()
    {scalar*_i=(scalar*)( p.b);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      { _stencil_val(s,0,0,0); }}}end_foreach_stencil();
  
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel reduction(+:sum)){
#line 164
foreach ()
    {scalar*_i=(scalar*)( p.b);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      sum += val(s,0,0,0);}}end_foreach();mpi_all_reduce_array(&sum,double,MPI_SUM,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
  
#line 167
s.sum = sum;
  s.nrelax = p.nrelax > 0 ? p.nrelax : 4;




  double resb;
  resb = s.resb = s.resa = p.residual (p.a, p.b, res, p.data);






  if (p.tolerance == 0.)
    p.tolerance = TOLERANCE;
  for (s.i = 0;
       s.i < NITERMAX && (s.i < NITERMIN || s.resa > p.tolerance);
       s.i++) {
    mg_cycle (p.a, res, da, p.relax, p.data,
       s.nrelax,
       p.minlevel,
       grid->maxdepth);
    s.resa = p.residual (p.a, p.b, res, p.data);
#line 199 "/u/basilisk/src/poisson.h"
    if (s.resa > p.tolerance) {
      if (resb/s.resa < 1.2 && s.nrelax < 100)
 s.nrelax++;
      else if (resb/s.resa > 10 && s.nrelax > 2)
 s.nrelax--;
    }







    resb = s.resa;
  }
  s.minlevel = p.minlevel;




  if (s.resa > p.tolerance) {
    scalar v = p.a[0];
    fprintf (ferr,
      "WARNING: convergence for %s not reached after %d iterations\n"
      "  res: %g sum: %g nrelax: %d\n", _attribute[v.i].name,
      s.i, s.resa, s.sum, s.nrelax), fflush (ferr);
  }




  if (!p.res)
    delete (res), pfree (res,__func__,__FILE__,0);
  delete (da), pfree (da,__func__,__FILE__,0);

  return s;
}
#line 258 "/u/basilisk/src/poisson.h"
struct Poisson {
  scalar a, b;
          vector alpha;
          scalar lambda;
  double tolerance;
  int nrelax, minlevel;
  scalar * res;

  double (* embed_flux) (Point, scalar, vector, double *);

};





static void relax (scalar * al, scalar * bl, int l, void * data)
{
  scalar a = al[0], b = bl[0];
  struct Poisson * p = (struct Poisson *) data;
          vector alpha = p->alpha;
          scalar lambda = p->lambda;
#line 296 "/u/basilisk/src/poisson.h"
  scalar c = a;






  if(!is_constant(lambda) && !is_constant(alpha.x)){{foreach_level_or_leaf (l) {
    double n = - sq(Delta)*val(b,0,0,0), d = - val(lambda,0,0,0)*sq(Delta);
     {
      n += val(alpha.x,1,0,0)*val(a,1,0,0) + val(alpha.x,0,0,0)*val(a,-1,0,0);
      d += val(alpha.x,1,0,0) + val(alpha.x,0,0,0);
    } 
#line 305
{
      n += val(alpha.y,0,1,0)*val(a,0,1,0) + val(alpha.y,0,0,0)*val(a,0,-1,0);
      d += val(alpha.y,0,1,0) + val(alpha.y,0,0,0);
    }

    if (p->embed_flux) {
      double c, e = p->embed_flux (point, a, alpha, &c);
      n -= c*sq(Delta);
      d += e*sq(Delta);
    }
    if (!d)
      val(c,0,0,0) = val(b,0,0,0) = 0.;
    else

      val(c,0,0,0) = n/d;
  }end_foreach_level_or_leaf();}}else if(is_constant(lambda) && !is_constant(alpha.x)){double _const_lambda=_constant[lambda.i-_NVARMAX];NOT_UNUSED(_const_lambda);






  {
#line 303
foreach_level_or_leaf (l) {
    double n = - sq(Delta)*val(b,0,0,0), d = - _const_lambda*sq(Delta);
     {
      n += val(alpha.x,1,0,0)*val(a,1,0,0) + val(alpha.x,0,0,0)*val(a,-1,0,0);
      d += val(alpha.x,1,0,0) + val(alpha.x,0,0,0);
    } 
#line 305
{
      n += val(alpha.y,0,1,0)*val(a,0,1,0) + val(alpha.y,0,0,0)*val(a,0,-1,0);
      d += val(alpha.y,0,1,0) + val(alpha.y,0,0,0);
    }

    if (p->embed_flux) {
      double c, e = p->embed_flux (point, a, alpha, &c);
      n -= c*sq(Delta);
      d += e*sq(Delta);
    }
    if (!d)
      val(c,0,0,0) = val(b,0,0,0) = 0.;
    else

      val(c,0,0,0) = n/d;
  }end_foreach_level_or_leaf();}}else if(!is_constant(lambda) && is_constant(alpha.x)){struct{double x,y;}_const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);






  {
#line 303
foreach_level_or_leaf (l) {
    double n = - sq(Delta)*val(b,0,0,0), d = - val(lambda,0,0,0)*sq(Delta);
     {
      n += _const_alpha.x*val(a,1,0,0) + _const_alpha.x*val(a,-1,0,0);
      d += _const_alpha.x + _const_alpha.x;
    } 
#line 305
{
      n += _const_alpha.y*val(a,0,1,0) + _const_alpha.y*val(a,0,-1,0);
      d += _const_alpha.y + _const_alpha.y;
    }

    if (p->embed_flux) {
      double c, e = p->embed_flux (point, a, alpha, &c);
      n -= c*sq(Delta);
      d += e*sq(Delta);
    }
    if (!d)
      val(c,0,0,0) = val(b,0,0,0) = 0.;
    else

      val(c,0,0,0) = n/d;
  }end_foreach_level_or_leaf();}}else {double _const_lambda=_constant[lambda.i-_NVARMAX];NOT_UNUSED(_const_lambda);struct{double x,y;}_const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);






  {
#line 303
foreach_level_or_leaf (l) {
    double n = - sq(Delta)*val(b,0,0,0), d = - _const_lambda*sq(Delta);
     {
      n += _const_alpha.x*val(a,1,0,0) + _const_alpha.x*val(a,-1,0,0);
      d += _const_alpha.x + _const_alpha.x;
    } 
#line 305
{
      n += _const_alpha.y*val(a,0,1,0) + _const_alpha.y*val(a,0,-1,0);
      d += _const_alpha.y + _const_alpha.y;
    }

    if (p->embed_flux) {
      double c, e = p->embed_flux (point, a, alpha, &c);
      n -= c*sq(Delta);
      d += e*sq(Delta);
    }
    if (!d)
      val(c,0,0,0) = val(b,0,0,0) = 0.;
    else

      val(c,0,0,0) = n/d;
  }end_foreach_level_or_leaf();}}
#line 338 "/u/basilisk/src/poisson.h"
}






static double residual (scalar * al, scalar * bl, scalar * resl, void * data)
{
  scalar a = al[0], b = bl[0], res = resl[0];
  struct Poisson * p = (struct Poisson *) data;
          vector alpha = p->alpha;
          scalar lambda = p->lambda;
  double maxres = 0.;


  vector  g=new_face_vector("g");
  foreach_face_stencil(){_stencil_is_face_x(){
    {_stencil_val_a(g.x,0,0,0); _stencil_val(alpha.x,0,0,0); _stencil_val(fs.x,0,0,0); _stencil_val(fs.x,0,0,0); _stencil_embed_face_gradient_x (point, a, 0);_stencil_val(a,0,0,0); _stencil_val(a,0 -1,0,0);           }}end__stencil_is_face_x()
#line 355
_stencil_is_face_y(){
    {_stencil_val_a(g.y,0,0,0); _stencil_val(alpha.y,0,0,0); _stencil_val(fs.y,0,0,0); _stencil_val(fs.y,0,0,0); _stencil_embed_face_gradient_y (point, a, 0);_stencil_val(a,0,0,0); _stencil_val(a,0,0 -1,0);           }}end__stencil_is_face_y()}end_foreach_face_stencil();
  
#line 355
if(!is_constant(alpha.x)){{foreach_face_generic(){is_face_x(){
    val(g.x,0,0,0) = val(alpha.x,0,0,0)*(_attribute[a.i].third && val(fs.x,0,0,0) < 1. && val(fs.x,0,0,0) > 0. ? embed_face_gradient_x (point, a, 0) : (val(a,0,0,0) - val(a,0 -1,0,0))/Delta);}end_is_face_x()
#line 355
is_face_y(){
    val(g.y,0,0,0) = val(alpha.y,0,0,0)*(_attribute[a.i].third && val(fs.y,0,0,0) < 1. && val(fs.y,0,0,0) > 0. ? embed_face_gradient_y (point, a, 0) : (val(a,0,0,0) - val(a,0,0 -1,0))/Delta);}end_is_face_y()}end_foreach_face_generic();}}else {struct{double x,y;}_const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);
  {
#line 355
foreach_face_generic(){is_face_x(){
    val(g.x,0,0,0) = _const_alpha.x*(_attribute[a.i].third && val(fs.x,0,0,0) < 1. && val(fs.x,0,0,0) > 0. ? embed_face_gradient_x (point, a, 0) : (val(a,0,0,0) - val(a,0 -1,0,0))/Delta);}end_is_face_x()
#line 355
is_face_y(){
    val(g.y,0,0,0) = _const_alpha.y*(_attribute[a.i].third && val(fs.y,0,0,0) < 1. && val(fs.y,0,0,0) > 0. ? embed_face_gradient_y (point, a, 0) : (val(a,0,0,0) - val(a,0,0 -1,0))/Delta);}end_is_face_y()}end_foreach_face_generic();}}
  foreach_stencil () {
    _stencil_val_a(res,0,0,0); _stencil_val(b,0,0,0); _stencil_val(lambda,0,0,0);_stencil_val(a,0,0,0);  
    
      {_stencil_val_r(res,0,0,0);_stencil_val(g.x,1,0,0); _stencil_val(g.x,0,0,0);   }
      
#line 360
{_stencil_val_r(res,0,0,0);_stencil_val(g.y,0,1,0); _stencil_val(g.y,0,0,0);   }

    if (p->embed_flux) {   
      default_stencil ( point,((scalar[]){ alpha.x,alpha.y, a,{-1}}) );
      _stencil_val_r(res,0,0,0);_stencil_val(a,0,0,0);    
    }

_stencil_val(res,0,0,0);
      {_stencil_val(res,0,0,0);   }

        
  
#line 369
}end_foreach_stencil();
  
#line 357
if(!is_constant(lambda)){
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel reduction(max:maxres)){
#line 357
foreach () {
    val(res,0,0,0) = val(b,0,0,0) - val(lambda,0,0,0)*val(a,0,0,0);
    
      val(res,0,0,0) -= (val(g.x,1,0,0) - val(g.x,0,0,0))/Delta;
      
#line 360
val(res,0,0,0) -= (val(g.y,0,1,0) - val(g.y,0,0,0))/Delta;

    if (p->embed_flux) {
      double c, e = p->embed_flux (point, a, alpha, &c);
      val(res,0,0,0) += c - e*val(a,0,0,0);
    }

    if (fabs (val(res,0,0,0)) > maxres)
      maxres = fabs (val(res,0,0,0));
  }end_foreach();mpi_all_reduce_array(&maxres,double,MPI_MAX,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 369
}else {double _const_lambda=_constant[lambda.i-_NVARMAX];NOT_UNUSED(_const_lambda);
  
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel reduction(max:maxres)){
#line 357
foreach () {
    val(res,0,0,0) = val(b,0,0,0) - _const_lambda*val(a,0,0,0);
    
      val(res,0,0,0) -= (val(g.x,1,0,0) - val(g.x,0,0,0))/Delta;
      
#line 360
val(res,0,0,0) -= (val(g.y,0,1,0) - val(g.y,0,0,0))/Delta;

    if (p->embed_flux) {
      double c, e = p->embed_flux (point, a, alpha, &c);
      val(res,0,0,0) += c - e*val(a,0,0,0);
    }

    if (fabs (val(res,0,0,0)) > maxres)
      maxres = fabs (val(res,0,0,0));
  }end_foreach();mpi_all_reduce_array(&maxres,double,MPI_MAX,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 369
}
#line 387 "/u/basilisk/src/poisson.h"
  {delete((scalar*)((vector[]){g,{{-1},{-1}}}));return maxres;}delete((scalar*)((vector[]){g,{{-1},{-1}}}));
}
#line 399 "/u/basilisk/src/poisson.h"
mgstats poisson (struct Poisson p)
{






  if (!p.alpha.x.i)
    p.alpha = unityf;
  if (!p.lambda.i)
    p.lambda = zeroc;




  vector alpha = p.alpha;
  scalar lambda = p.lambda;
  restriction (((scalar[]){alpha.x,alpha.y,lambda,{-1}}));





  double defaultol = TOLERANCE;
  if (p.tolerance)
    TOLERANCE = p.tolerance;

  scalar a = p.a, b = p.b;

  if (!p.embed_flux && _attribute[a.i].boundary[embed] != symmetry)
    p.embed_flux = embed_flux;

  mgstats s = mg_solve ((struct MGSolve){((scalar[]){a,{-1}}), ((scalar[]){b,{-1}}), residual, relax,
   &p, p.nrelax, p.res, .minlevel = max(1, p.minlevel)});




  if (p.tolerance)
    TOLERANCE = defaultol;

  return s;
}
#line 461 "/u/basilisk/src/poisson.h"
struct Project {
  vector uf;
  scalar p;
  vector alpha;
  double dt;
  int nrelax;
};

     
mgstats project (struct Project q)
{tracing("project","/u/basilisk/src/poisson.h",0);
  vector uf = q.uf;
  scalar p = q.p;
          vector alpha = q.alpha.x.i ? q.alpha : unityf;
  double dt = q.dt ? q.dt : 1.;
  int nrelax = q.nrelax ? q.nrelax : 4;






  scalar  div=new_scalar("div");
  foreach_stencil() {
    _stencil_val_a(div,0,0,0);  
    
      {_stencil_val_r(div,0,0,0); _stencil_val(uf.x,1,0,0); _stencil_val(uf.x,0,0,0);  }
      
#line 487
{_stencil_val_r(div,0,0,0); _stencil_val(uf.y,0,1,0); _stencil_val(uf.y,0,0,0);  }
    _stencil_val_r(div,0,0,0);  
  }end_foreach_stencil();
  {
#line 484
foreach() {
    val(div,0,0,0) = 0.;
    
      val(div,0,0,0) += val(uf.x,1,0,0) - val(uf.x,0,0,0);
      
#line 487
val(div,0,0,0) += val(uf.y,0,1,0) - val(uf.y,0,0,0);
    val(div,0,0,0) /= dt*Delta;
  }end_foreach();}
#line 500 "/u/basilisk/src/poisson.h"
  mgstats mgp = poisson ((struct Poisson){p, div, alpha,
    .tolerance = TOLERANCE/sq(dt), .nrelax = nrelax});




  foreach_face_stencil(){_stencil_is_face_x(){
    {_stencil_val_r(uf.x,0,0,0);_stencil_val(alpha.x,0,0,0); _stencil_val(fs.x,0,0,0); _stencil_val(fs.x,0,0,0); _stencil_embed_face_gradient_x (point, p, 0);_stencil_val(p,0,0,0); _stencil_val(p,0 -1,0,0);            }}end__stencil_is_face_x()
#line 506
_stencil_is_face_y(){
    {_stencil_val_r(uf.y,0,0,0);_stencil_val(alpha.y,0,0,0); _stencil_val(fs.y,0,0,0); _stencil_val(fs.y,0,0,0); _stencil_embed_face_gradient_y (point, p, 0);_stencil_val(p,0,0,0); _stencil_val(p,0,0 -1,0);            }}end__stencil_is_face_y()}end_foreach_face_stencil();




  
#line 506
if(!is_constant(alpha.x)){{foreach_face_generic(){is_face_x(){
    val(uf.x,0,0,0) -= dt*val(alpha.x,0,0,0)*(_attribute[p.i].third && val(fs.x,0,0,0) < 1. && val(fs.x,0,0,0) > 0. ? embed_face_gradient_x (point, p, 0) : (val(p,0,0,0) - val(p,0 -1,0,0))/Delta);}end_is_face_x()
#line 506
is_face_y(){
    val(uf.y,0,0,0) -= dt*val(alpha.y,0,0,0)*(_attribute[p.i].third && val(fs.y,0,0,0) < 1. && val(fs.y,0,0,0) > 0. ? embed_face_gradient_y (point, p, 0) : (val(p,0,0,0) - val(p,0,0 -1,0))/Delta);}end_is_face_y()}end_foreach_face_generic();}}else {struct{double x,y;}_const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);




  {
#line 506
foreach_face_generic(){is_face_x(){
    val(uf.x,0,0,0) -= dt*_const_alpha.x*(_attribute[p.i].third && val(fs.x,0,0,0) < 1. && val(fs.x,0,0,0) > 0. ? embed_face_gradient_x (point, p, 0) : (val(p,0,0,0) - val(p,0 -1,0,0))/Delta);}end_is_face_x()
#line 506
is_face_y(){
    val(uf.y,0,0,0) -= dt*_const_alpha.y*(_attribute[p.i].third && val(fs.y,0,0,0) < 1. && val(fs.y,0,0,0) > 0. ? embed_face_gradient_y (point, p, 0) : (val(p,0,0,0) - val(p,0,0 -1,0))/Delta);}end_is_face_y()}end_foreach_face_generic();}}

  {delete((scalar*)((scalar[]){div,{-1}}));{end_tracing("project","/u/basilisk/src/poisson.h",0);return mgp;}}delete((scalar*)((scalar[]){div,{-1}}));
end_tracing("project","/u/basilisk/src/poisson.h",0);}
#line 2 "/u/basilisk/src/viscosity-embed.h"

struct Viscosity {
  vector u;
  vector mu;
  scalar rho;
  double dt;
  int nrelax;
  scalar * res;
  double (* embed_flux) (Point, scalar, vector, double *);
};
#line 30 "/u/basilisk/src/viscosity-embed.h"
static void relax_diffusion (scalar * a, scalar * b, int l, void * data)
{
  struct Viscosity * p = (struct Viscosity *) data;
          vector mu = p->mu;
          scalar rho = p->rho;
  double dt = p->dt;
  vector u = (*((vector *)&(a[0]))), r = (*((vector *)&(b[0])));

  double (* embed_flux) (Point, scalar, vector, double *) = p->embed_flux;
  if(!is_constant(mu.x) && !is_constant(rho)){{foreach_level_or_leaf (l) {
    double avgmu = 0.;
    
      avgmu += val(mu.x,0,0,0) + val(mu.x,1,0,0);
      
#line 42
avgmu += val(mu.y,0,0,0) + val(mu.y,0,1,0);
    avgmu = dt*avgmu + 1e-30;
     {
      double c = 0.;
      double d = embed_flux ? embed_flux (point, u.x, mu, &c) : 0.;
      scalar s = u.x;
      double a = 0.;
      
 a += val(mu.x,1,0,0)*val(s,1,0,0) + val(mu.x,0,0,0)*val(s,-1,0,0);
 
#line 50
a += val(mu.y,0,1,0)*val(s,0,1,0) + val(mu.y,0,0,0)*val(s,0,-1,0);
      val(u.x,0,0,0) = (dt*a + (val(r.x,0,0,0) - dt*c)*sq(Delta))/
 (sq(Delta)*(val(rho,0,0,0)*((coord){1.,1.}).x + dt*d) + avgmu);
    } 
#line 44
{
      double c = 0.;
      double d = embed_flux ? embed_flux (point, u.y, mu, &c) : 0.;
      scalar s = u.y;
      double a = 0.;
      
 a += val(mu.y,0,1,0)*val(s,0,1,0) + val(mu.y,0,0,0)*val(s,0,-1,0);
 
#line 50
a += val(mu.x,1,0,0)*val(s,1,0,0) + val(mu.x,0,0,0)*val(s,-1,0,0);
      val(u.y,0,0,0) = (dt*a + (val(r.y,0,0,0) - dt*c)*sq(Delta))/
 (sq(Delta)*(val(rho,0,0,0)*((coord){1.,1.}).y + dt*d) + avgmu);
    }
  }end_foreach_level_or_leaf();}}else if(is_constant(mu.x) && !is_constant(rho)){struct{double x,y;}_const_mu={_constant[mu.x.i-_NVARMAX],_constant[mu.y.i-_NVARMAX]};NOT_UNUSED(_const_mu);
  {
#line 39
foreach_level_or_leaf (l) {
    double avgmu = 0.;
    
      avgmu += _const_mu.x + _const_mu.x;
      
#line 42
avgmu += _const_mu.y + _const_mu.y;
    avgmu = dt*avgmu + 1e-30;
     {
      double c = 0.;
      double d = embed_flux ? embed_flux (point, u.x, mu, &c) : 0.;
      scalar s = u.x;
      double a = 0.;
      
 a += _const_mu.x*val(s,1,0,0) + _const_mu.x*val(s,-1,0,0);
 
#line 50
a += _const_mu.y*val(s,0,1,0) + _const_mu.y*val(s,0,-1,0);
      val(u.x,0,0,0) = (dt*a + (val(r.x,0,0,0) - dt*c)*sq(Delta))/
 (sq(Delta)*(val(rho,0,0,0)*((coord){1.,1.}).x + dt*d) + avgmu);
    } 
#line 44
{
      double c = 0.;
      double d = embed_flux ? embed_flux (point, u.y, mu, &c) : 0.;
      scalar s = u.y;
      double a = 0.;
      
 a += _const_mu.y*val(s,0,1,0) + _const_mu.y*val(s,0,-1,0);
 
#line 50
a += _const_mu.x*val(s,1,0,0) + _const_mu.x*val(s,-1,0,0);
      val(u.y,0,0,0) = (dt*a + (val(r.y,0,0,0) - dt*c)*sq(Delta))/
 (sq(Delta)*(val(rho,0,0,0)*((coord){1.,1.}).y + dt*d) + avgmu);
    }
  }end_foreach_level_or_leaf();}}else if(!is_constant(mu.x) && is_constant(rho)){double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);
  {
#line 39
foreach_level_or_leaf (l) {
    double avgmu = 0.;
    
      avgmu += val(mu.x,0,0,0) + val(mu.x,1,0,0);
      
#line 42
avgmu += val(mu.y,0,0,0) + val(mu.y,0,1,0);
    avgmu = dt*avgmu + 1e-30;
     {
      double c = 0.;
      double d = embed_flux ? embed_flux (point, u.x, mu, &c) : 0.;
      scalar s = u.x;
      double a = 0.;
      
 a += val(mu.x,1,0,0)*val(s,1,0,0) + val(mu.x,0,0,0)*val(s,-1,0,0);
 
#line 50
a += val(mu.y,0,1,0)*val(s,0,1,0) + val(mu.y,0,0,0)*val(s,0,-1,0);
      val(u.x,0,0,0) = (dt*a + (val(r.x,0,0,0) - dt*c)*sq(Delta))/
 (sq(Delta)*(_const_rho*((coord){1.,1.}).x + dt*d) + avgmu);
    } 
#line 44
{
      double c = 0.;
      double d = embed_flux ? embed_flux (point, u.y, mu, &c) : 0.;
      scalar s = u.y;
      double a = 0.;
      
 a += val(mu.y,0,1,0)*val(s,0,1,0) + val(mu.y,0,0,0)*val(s,0,-1,0);
 
#line 50
a += val(mu.x,1,0,0)*val(s,1,0,0) + val(mu.x,0,0,0)*val(s,-1,0,0);
      val(u.y,0,0,0) = (dt*a + (val(r.y,0,0,0) - dt*c)*sq(Delta))/
 (sq(Delta)*(_const_rho*((coord){1.,1.}).y + dt*d) + avgmu);
    }
  }end_foreach_level_or_leaf();}}else {struct{double x,y;}_const_mu={_constant[mu.x.i-_NVARMAX],_constant[mu.y.i-_NVARMAX]};NOT_UNUSED(_const_mu);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);
  {
#line 39
foreach_level_or_leaf (l) {
    double avgmu = 0.;
    
      avgmu += _const_mu.x + _const_mu.x;
      
#line 42
avgmu += _const_mu.y + _const_mu.y;
    avgmu = dt*avgmu + 1e-30;
     {
      double c = 0.;
      double d = embed_flux ? embed_flux (point, u.x, mu, &c) : 0.;
      scalar s = u.x;
      double a = 0.;
      
 a += _const_mu.x*val(s,1,0,0) + _const_mu.x*val(s,-1,0,0);
 
#line 50
a += _const_mu.y*val(s,0,1,0) + _const_mu.y*val(s,0,-1,0);
      val(u.x,0,0,0) = (dt*a + (val(r.x,0,0,0) - dt*c)*sq(Delta))/
 (sq(Delta)*(_const_rho*((coord){1.,1.}).x + dt*d) + avgmu);
    } 
#line 44
{
      double c = 0.;
      double d = embed_flux ? embed_flux (point, u.y, mu, &c) : 0.;
      scalar s = u.y;
      double a = 0.;
      
 a += _const_mu.y*val(s,0,1,0) + _const_mu.y*val(s,0,-1,0);
 
#line 50
a += _const_mu.x*val(s,1,0,0) + _const_mu.x*val(s,-1,0,0);
      val(u.y,0,0,0) = (dt*a + (val(r.y,0,0,0) - dt*c)*sq(Delta))/
 (sq(Delta)*(_const_rho*((coord){1.,1.}).y + dt*d) + avgmu);
    }
  }end_foreach_level_or_leaf();}}
#line 66 "/u/basilisk/src/viscosity-embed.h"
}

static double residual_diffusion (scalar * a, scalar * b, scalar * resl,
      void * data)
{
  struct Viscosity * p = (struct Viscosity *) data;
          vector mu = p->mu;
          scalar rho = p->rho;
  double dt = p->dt;
  double (* embed_flux) (Point, scalar, vector, double *) = p->embed_flux;
  vector u = (*((vector *)&(a[0]))), r = (*((vector *)&(b[0]))), res = (*((vector *)&(resl[0])));
  double maxres = 0.;


   {
    scalar s = u.x;
    vector  g=new_face_vector("g");
    foreach_face_stencil(){_stencil_is_face_x(){
      {_stencil_val_a(g.x,0,0,0); _stencil_val(mu.x,0,0,0); _stencil_val(fs.x,0,0,0); _stencil_val(fs.x,0,0,0); _stencil_embed_face_gradient_x (point, s, 0);_stencil_val(s,0,0,0); _stencil_val(s,0 -1,0,0);           }}end__stencil_is_face_x()
#line 83
_stencil_is_face_y(){
      {_stencil_val_a(g.y,0,0,0); _stencil_val(mu.y,0,0,0); _stencil_val(fs.y,0,0,0); _stencil_val(fs.y,0,0,0); _stencil_embed_face_gradient_y (point, s, 0);_stencil_val(s,0,0,0); _stencil_val(s,0,0 -1,0);           }}end__stencil_is_face_y()}end_foreach_face_stencil();
    
#line 83
if(!is_constant(mu.x)){{foreach_face_generic(){is_face_x(){
      val(g.x,0,0,0) = val(mu.x,0,0,0)*(_attribute[s.i].third && val(fs.x,0,0,0) < 1. && val(fs.x,0,0,0) > 0. ? embed_face_gradient_x (point, s, 0) : (val(s,0,0,0) - val(s,0 -1,0,0))/Delta);}end_is_face_x()
#line 83
is_face_y(){
      val(g.y,0,0,0) = val(mu.y,0,0,0)*(_attribute[s.i].third && val(fs.y,0,0,0) < 1. && val(fs.y,0,0,0) > 0. ? embed_face_gradient_y (point, s, 0) : (val(s,0,0,0) - val(s,0,0 -1,0))/Delta);}end_is_face_y()}end_foreach_face_generic();}}else {struct{double x,y;}_const_mu={_constant[mu.x.i-_NVARMAX],_constant[mu.y.i-_NVARMAX]};NOT_UNUSED(_const_mu);
    {
#line 83
foreach_face_generic(){is_face_x(){
      val(g.x,0,0,0) = _const_mu.x*(_attribute[s.i].third && val(fs.x,0,0,0) < 1. && val(fs.x,0,0,0) > 0. ? embed_face_gradient_x (point, s, 0) : (val(s,0,0,0) - val(s,0 -1,0,0))/Delta);}end_is_face_x()
#line 83
is_face_y(){
      val(g.y,0,0,0) = _const_mu.y*(_attribute[s.i].third && val(fs.y,0,0,0) < 1. && val(fs.y,0,0,0) > 0. ? embed_face_gradient_y (point, s, 0) : (val(s,0,0,0) - val(s,0,0 -1,0))/Delta);}end_is_face_y()}end_foreach_face_generic();}}
    foreach_stencil () {   
      
      
 { _stencil_val(g.x,0,0,0); _stencil_val(g.x,1,0,0);  }
 
#line 88
{ _stencil_val(g.y,0,0,0); _stencil_val(g.y,0,1,0);  }
      _stencil_val_a(res.x,0,0,0); _stencil_val(r.x,0,0,0); _stencil_val(rho,0,0,0);_stencil_val(u.x,0,0,0);    
      if (embed_flux) {   
 default_stencil ( point,((scalar[]){ mu.x,mu.y, u.x,{-1}}) );
 _stencil_val_r(res.x,0,0,0);_stencil_val(u.x,0,0,0);    
      }
_stencil_val(res.x,0,0,0);
 {_stencil_val(res.x,0,0,0);   }
          
    
#line 96
}end_foreach_stencil();
    
#line 85
if(!is_constant(rho)){
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel reduction(max:maxres)){
#line 85
foreach () {
      double a = 0.;
      
 a += val(g.x,0,0,0) - val(g.x,1,0,0);
 
#line 88
a += val(g.y,0,0,0) - val(g.y,0,1,0);
      val(res.x,0,0,0) = val(r.x,0,0,0) - val(rho,0,0,0)*((coord){1.,1.}).x*val(u.x,0,0,0) - dt*a/Delta;
      if (embed_flux) {
 double c, d = embed_flux (point, u.x, mu, &c);
 val(res.x,0,0,0) -= dt*(c + d*val(u.x,0,0,0));
      }
      if (fabs (val(res.x,0,0,0)) > maxres)
 maxres = fabs (val(res.x,0,0,0));
    }end_foreach();mpi_all_reduce_array(&maxres,double,MPI_MAX,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 96
}else {double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);
    
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel reduction(max:maxres)){
#line 85
foreach () {
      double a = 0.;
      
 a += val(g.x,0,0,0) - val(g.x,1,0,0);
 
#line 88
a += val(g.y,0,0,0) - val(g.y,0,1,0);
      val(res.x,0,0,0) = val(r.x,0,0,0) - _const_rho*((coord){1.,1.}).x*val(u.x,0,0,0) - dt*a/Delta;
      if (embed_flux) {
 double c, d = embed_flux (point, u.x, mu, &c);
 val(res.x,0,0,0) -= dt*(c + d*val(u.x,0,0,0));
      }
      if (fabs (val(res.x,0,0,0)) > maxres)
 maxres = fabs (val(res.x,0,0,0));
    }end_foreach();mpi_all_reduce_array(&maxres,double,MPI_MAX,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 96
}delete((scalar*)((vector[]){g,{{-1},{-1}}}));
  } 
#line 80
{
    scalar s = u.y;
    vector  g=new_face_vector("g");
    foreach_face_stencil(){_stencil_is_face_y(){
      {_stencil_val_a(g.y,0,0,0); _stencil_val(mu.y,0,0,0); _stencil_val(fs.y,0,0,0); _stencil_val(fs.y,0,0,0); _stencil_embed_face_gradient_y (point, s, 0);_stencil_val(s,0,0,0); _stencil_val(s,0,0 -1,0);           }}end__stencil_is_face_y()
#line 83
_stencil_is_face_x(){
      {_stencil_val_a(g.x,0,0,0); _stencil_val(mu.x,0,0,0); _stencil_val(fs.x,0,0,0); _stencil_val(fs.x,0,0,0); _stencil_embed_face_gradient_x (point, s, 0);_stencil_val(s,0,0,0); _stencil_val(s,0 -1,0,0);           }}end__stencil_is_face_x()}end_foreach_face_stencil();
    
#line 83
if(!is_constant(mu.y)){{foreach_face_generic(){is_face_y(){
      val(g.y,0,0,0) = val(mu.y,0,0,0)*(_attribute[s.i].third && val(fs.y,0,0,0) < 1. && val(fs.y,0,0,0) > 0. ? embed_face_gradient_y (point, s, 0) : (val(s,0,0,0) - val(s,0,0 -1,0))/Delta);}end_is_face_y()
#line 83
is_face_x(){
      val(g.x,0,0,0) = val(mu.x,0,0,0)*(_attribute[s.i].third && val(fs.x,0,0,0) < 1. && val(fs.x,0,0,0) > 0. ? embed_face_gradient_x (point, s, 0) : (val(s,0,0,0) - val(s,0 -1,0,0))/Delta);}end_is_face_x()}end_foreach_face_generic();}}else {struct{double x,y;}_const_mu={_constant[mu.y.i-_NVARMAX],_constant[mu.x.i-_NVARMAX]};NOT_UNUSED(_const_mu);
    {
#line 83
foreach_face_generic(){is_face_y(){
      val(g.y,0,0,0) = _const_mu.y*(_attribute[s.i].third && val(fs.y,0,0,0) < 1. && val(fs.y,0,0,0) > 0. ? embed_face_gradient_y (point, s, 0) : (val(s,0,0,0) - val(s,0,0 -1,0))/Delta);}end_is_face_y()
#line 83
is_face_x(){
      val(g.x,0,0,0) = _const_mu.x*(_attribute[s.i].third && val(fs.x,0,0,0) < 1. && val(fs.x,0,0,0) > 0. ? embed_face_gradient_x (point, s, 0) : (val(s,0,0,0) - val(s,0 -1,0,0))/Delta);}end_is_face_x()}end_foreach_face_generic();}}
    foreach_stencil () {   
      
      
 { _stencil_val(g.y,0,0,0); _stencil_val(g.y,0,1,0);  }
 
#line 88
{ _stencil_val(g.x,0,0,0); _stencil_val(g.x,1,0,0);  }
      _stencil_val_a(res.y,0,0,0); _stencil_val(r.y,0,0,0); _stencil_val(rho,0,0,0);_stencil_val(u.y,0,0,0);    
      if (embed_flux) {   
 default_stencil ( point,((scalar[]){ mu.x,mu.y, u.y,{-1}}) );
 _stencil_val_r(res.y,0,0,0);_stencil_val(u.y,0,0,0);    
      }
_stencil_val(res.y,0,0,0);
 {_stencil_val(res.y,0,0,0);   }
          
    
#line 96
}end_foreach_stencil();
    
#line 85
if(!is_constant(rho)){
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel reduction(max:maxres)){
#line 85
foreach () {
      double a = 0.;
      
 a += val(g.y,0,0,0) - val(g.y,0,1,0);
 
#line 88
a += val(g.x,0,0,0) - val(g.x,1,0,0);
      val(res.y,0,0,0) = val(r.y,0,0,0) - val(rho,0,0,0)*((coord){1.,1.}).y*val(u.y,0,0,0) - dt*a/Delta;
      if (embed_flux) {
 double c, d = embed_flux (point, u.y, mu, &c);
 val(res.y,0,0,0) -= dt*(c + d*val(u.y,0,0,0));
      }
      if (fabs (val(res.y,0,0,0)) > maxres)
 maxres = fabs (val(res.y,0,0,0));
    }end_foreach();mpi_all_reduce_array(&maxres,double,MPI_MAX,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 96
}else {double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);
    
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel reduction(max:maxres)){
#line 85
foreach () {
      double a = 0.;
      
 a += val(g.y,0,0,0) - val(g.y,0,1,0);
 
#line 88
a += val(g.x,0,0,0) - val(g.x,1,0,0);
      val(res.y,0,0,0) = val(r.y,0,0,0) - _const_rho*((coord){1.,1.}).y*val(u.y,0,0,0) - dt*a/Delta;
      if (embed_flux) {
 double c, d = embed_flux (point, u.y, mu, &c);
 val(res.y,0,0,0) -= dt*(c + d*val(u.y,0,0,0));
      }
      if (fabs (val(res.y,0,0,0)) > maxres)
 maxres = fabs (val(res.y,0,0,0));
    }end_foreach();mpi_all_reduce_array(&maxres,double,MPI_MAX,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 96
}delete((scalar*)((vector[]){g,{{-1},{-1}}}));
  }
#line 115 "/u/basilisk/src/viscosity-embed.h"
  return maxres;
}



double TOLERANCE_MU = 0.;

     
mgstats viscosity (struct Viscosity p)
{tracing("viscosity","/u/basilisk/src/viscosity-embed.h",0);
  vector u = p.u,  r=new_vector("r");
  scalar rho = p.rho;
  foreach_stencil()
    {
      {_stencil_val_a(r.x,0,0,0); _stencil_val(rho,0,0,0);_stencil_val(u.x,0,0,0); }
      
#line 129
{_stencil_val_a(r.y,0,0,0); _stencil_val(rho,0,0,0);_stencil_val(u.y,0,0,0); }}end_foreach_stencil();
  {
#line 127
foreach()
    {
      val(r.x,0,0,0) = val(rho,0,0,0)*val(u.x,0,0,0);
      
#line 129
val(r.y,0,0,0) = val(rho,0,0,0)*val(u.y,0,0,0);}end_foreach();}

  vector mu = p.mu;
  restriction (((scalar[]){mu.x,mu.y, rho,{-1}}));

  p.embed_flux = _attribute[u.x.i].boundary[embed] != antisymmetry ? embed_flux : NULL;
  { mgstats _ret= mg_solve ((struct MGSolve){(scalar *)((vector[]){u,{{-1},{-1}}}), (scalar *)((vector[]){r,{{-1},{-1}}}),
     residual_diffusion, relax_diffusion, &p, p.nrelax, p.res,
     .minlevel = 1,

     .tolerance = TOLERANCE_MU});delete((scalar*)((vector[]){r,{{-1},{-1}}}));{end_tracing("viscosity","/u/basilisk/src/viscosity-embed.h",0);return _ret;}}delete((scalar*)((vector[]){r,{{-1},{-1}}}));
end_tracing("viscosity","/u/basilisk/src/viscosity-embed.h",0);}
#line 32 "/u/basilisk/src/navier-stokes/centered.h"
#line 44 "/u/basilisk/src/navier-stokes/centered.h"
scalar  p={3};
vector  u={{4},{5}},  g={{6},{7}};
scalar  pf={8};
vector  uf={{9},{10}};
#line 70 "/u/basilisk/src/navier-stokes/centered.h"
        vector mu = {{_NVARMAX+0},{_NVARMAX+1}}, a = {{_NVARMAX+0},{_NVARMAX+1}}, alpha = {{_NVARMAX+2},{_NVARMAX+3}};
        scalar rho = {_NVARMAX+4};
mgstats mgp, mgpf, mgu;
bool stokes = false;
#line 91
static double _boundary0(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;if(!is_constant(alpha.x) && !is_constant(a.x) && !is_constant(fm.x) && !is_constant(rho) && !is_constant(cm)){{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann ((val(alpha.x,1,0,0) ? val(a.x,1,0,0)*val(fm.x,1,0,0)/val(alpha.x,1,0,0) : val(a.x,1,0,0)*val(rho,0,0,0)/(val(cm,0,0,0) + 1e-30)));}}}else if(is_constant(alpha.x) && !is_constant(a.x) && !is_constant(fm.x) && !is_constant(rho) && !is_constant(cm)){struct{double x,y;}_const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann ((_const_alpha.x ? val(a.x,1,0,0)*val(fm.x,1,0,0)/_const_alpha.x : val(a.x,1,0,0)*val(rho,0,0,0)/(val(cm,0,0,0) + 1e-30)));}}}else if(!is_constant(alpha.x) && is_constant(a.x) && !is_constant(fm.x) && !is_constant(rho) && !is_constant(cm)){struct{double x,y;}_const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann ((val(alpha.x,1,0,0) ? _const_a.x*val(fm.x,1,0,0)/val(alpha.x,1,0,0) : _const_a.x*val(rho,0,0,0)/(val(cm,0,0,0) + 1e-30)));}}}else if(is_constant(alpha.x) && is_constant(a.x) && !is_constant(fm.x) && !is_constant(rho) && !is_constant(cm)){struct{double x,y;}_const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);struct{double x,y;}_const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann ((_const_alpha.x ? _const_a.x*val(fm.x,1,0,0)/_const_alpha.x : _const_a.x*val(rho,0,0,0)/(val(cm,0,0,0) + 1e-30)));}}}else if(!is_constant(alpha.x) && !is_constant(a.x) && is_constant(fm.x) && !is_constant(rho) && !is_constant(cm)){struct{double x,y;}_const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann ((val(alpha.x,1,0,0) ? val(a.x,1,0,0)*_const_fm.x/val(alpha.x,1,0,0) : val(a.x,1,0,0)*val(rho,0,0,0)/(val(cm,0,0,0) + 1e-30)));}}}else if(is_constant(alpha.x) && !is_constant(a.x) && is_constant(fm.x) && !is_constant(rho) && !is_constant(cm)){struct{double x,y;}_const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);struct{double x,y;}_const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann ((_const_alpha.x ? val(a.x,1,0,0)*_const_fm.x/_const_alpha.x : val(a.x,1,0,0)*val(rho,0,0,0)/(val(cm,0,0,0) + 1e-30)));}}}else if(!is_constant(alpha.x) && is_constant(a.x) && is_constant(fm.x) && !is_constant(rho) && !is_constant(cm)){struct{double x,y;}_const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);struct{double x,y;}_const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann ((val(alpha.x,1,0,0) ? _const_a.x*_const_fm.x/val(alpha.x,1,0,0) : _const_a.x*val(rho,0,0,0)/(val(cm,0,0,0) + 1e-30)));}}}else if(is_constant(alpha.x) && is_constant(a.x) && is_constant(fm.x) && !is_constant(rho) && !is_constant(cm)){struct{double x,y;}_const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);struct{double x,y;}_const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);struct{double x,y;}_const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann ((_const_alpha.x ? _const_a.x*_const_fm.x/_const_alpha.x : _const_a.x*val(rho,0,0,0)/(val(cm,0,0,0) + 1e-30)));}}}else if(!is_constant(alpha.x) && !is_constant(a.x) && !is_constant(fm.x) && is_constant(rho) && !is_constant(cm)){double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann ((val(alpha.x,1,0,0) ? val(a.x,1,0,0)*val(fm.x,1,0,0)/val(alpha.x,1,0,0) : val(a.x,1,0,0)*_const_rho/(val(cm,0,0,0) + 1e-30)));}}}else if(is_constant(alpha.x) && !is_constant(a.x) && !is_constant(fm.x) && is_constant(rho) && !is_constant(cm)){struct{double x,y;}_const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann ((_const_alpha.x ? val(a.x,1,0,0)*val(fm.x,1,0,0)/_const_alpha.x : val(a.x,1,0,0)*_const_rho/(val(cm,0,0,0) + 1e-30)));}}}else if(!is_constant(alpha.x) && is_constant(a.x) && !is_constant(fm.x) && is_constant(rho) && !is_constant(cm)){struct{double x,y;}_const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann ((val(alpha.x,1,0,0) ? _const_a.x*val(fm.x,1,0,0)/val(alpha.x,1,0,0) : _const_a.x*_const_rho/(val(cm,0,0,0) + 1e-30)));}}}else if(is_constant(alpha.x) && is_constant(a.x) && !is_constant(fm.x) && is_constant(rho) && !is_constant(cm)){struct{double x,y;}_const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);struct{double x,y;}_const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann ((_const_alpha.x ? _const_a.x*val(fm.x,1,0,0)/_const_alpha.x : _const_a.x*_const_rho/(val(cm,0,0,0) + 1e-30)));}}}else if(!is_constant(alpha.x) && !is_constant(a.x) && is_constant(fm.x) && is_constant(rho) && !is_constant(cm)){struct{double x,y;}_const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann ((val(alpha.x,1,0,0) ? val(a.x,1,0,0)*_const_fm.x/val(alpha.x,1,0,0) : val(a.x,1,0,0)*_const_rho/(val(cm,0,0,0) + 1e-30)));}}}else if(is_constant(alpha.x) && !is_constant(a.x) && is_constant(fm.x) && is_constant(rho) && !is_constant(cm)){struct{double x,y;}_const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);struct{double x,y;}_const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann ((_const_alpha.x ? val(a.x,1,0,0)*_const_fm.x/_const_alpha.x : val(a.x,1,0,0)*_const_rho/(val(cm,0,0,0) + 1e-30)));}}}else if(!is_constant(alpha.x) && is_constant(a.x) && is_constant(fm.x) && is_constant(rho) && !is_constant(cm)){struct{double x,y;}_const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);struct{double x,y;}_const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann ((val(alpha.x,1,0,0) ? _const_a.x*_const_fm.x/val(alpha.x,1,0,0) : _const_a.x*_const_rho/(val(cm,0,0,0) + 1e-30)));}}}else if(is_constant(alpha.x) && is_constant(a.x) && is_constant(fm.x) && is_constant(rho) && !is_constant(cm)){struct{double x,y;}_const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);struct{double x,y;}_const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);struct{double x,y;}_const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann ((_const_alpha.x ? _const_a.x*_const_fm.x/_const_alpha.x : _const_a.x*_const_rho/(val(cm,0,0,0) + 1e-30)));}}}else if(!is_constant(alpha.x) && !is_constant(a.x) && !is_constant(fm.x) && !is_constant(rho) && is_constant(cm)){double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann ((val(alpha.x,1,0,0) ? val(a.x,1,0,0)*val(fm.x,1,0,0)/val(alpha.x,1,0,0) : val(a.x,1,0,0)*val(rho,0,0,0)/(_const_cm + 1e-30)));}}}else if(is_constant(alpha.x) && !is_constant(a.x) && !is_constant(fm.x) && !is_constant(rho) && is_constant(cm)){struct{double x,y;}_const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann ((_const_alpha.x ? val(a.x,1,0,0)*val(fm.x,1,0,0)/_const_alpha.x : val(a.x,1,0,0)*val(rho,0,0,0)/(_const_cm + 1e-30)));}}}else if(!is_constant(alpha.x) && is_constant(a.x) && !is_constant(fm.x) && !is_constant(rho) && is_constant(cm)){struct{double x,y;}_const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann ((val(alpha.x,1,0,0) ? _const_a.x*val(fm.x,1,0,0)/val(alpha.x,1,0,0) : _const_a.x*val(rho,0,0,0)/(_const_cm + 1e-30)));}}}else if(is_constant(alpha.x) && is_constant(a.x) && !is_constant(fm.x) && !is_constant(rho) && is_constant(cm)){struct{double x,y;}_const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);struct{double x,y;}_const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann ((_const_alpha.x ? _const_a.x*val(fm.x,1,0,0)/_const_alpha.x : _const_a.x*val(rho,0,0,0)/(_const_cm + 1e-30)));}}}else if(!is_constant(alpha.x) && !is_constant(a.x) && is_constant(fm.x) && !is_constant(rho) && is_constant(cm)){struct{double x,y;}_const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann ((val(alpha.x,1,0,0) ? val(a.x,1,0,0)*_const_fm.x/val(alpha.x,1,0,0) : val(a.x,1,0,0)*val(rho,0,0,0)/(_const_cm + 1e-30)));}}}else if(is_constant(alpha.x) && !is_constant(a.x) && is_constant(fm.x) && !is_constant(rho) && is_constant(cm)){struct{double x,y;}_const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);struct{double x,y;}_const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann ((_const_alpha.x ? val(a.x,1,0,0)*_const_fm.x/_const_alpha.x : val(a.x,1,0,0)*val(rho,0,0,0)/(_const_cm + 1e-30)));}}}else if(!is_constant(alpha.x) && is_constant(a.x) && is_constant(fm.x) && !is_constant(rho) && is_constant(cm)){struct{double x,y;}_const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);struct{double x,y;}_const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann ((val(alpha.x,1,0,0) ? _const_a.x*_const_fm.x/val(alpha.x,1,0,0) : _const_a.x*val(rho,0,0,0)/(_const_cm + 1e-30)));}}}else if(is_constant(alpha.x) && is_constant(a.x) && is_constant(fm.x) && !is_constant(rho) && is_constant(cm)){struct{double x,y;}_const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);struct{double x,y;}_const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);struct{double x,y;}_const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann ((_const_alpha.x ? _const_a.x*_const_fm.x/_const_alpha.x : _const_a.x*val(rho,0,0,0)/(_const_cm + 1e-30)));}}}else if(!is_constant(alpha.x) && !is_constant(a.x) && !is_constant(fm.x) && is_constant(rho) && is_constant(cm)){double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann ((val(alpha.x,1,0,0) ? val(a.x,1,0,0)*val(fm.x,1,0,0)/val(alpha.x,1,0,0) : val(a.x,1,0,0)*_const_rho/(_const_cm + 1e-30)));}}}else if(is_constant(alpha.x) && !is_constant(a.x) && !is_constant(fm.x) && is_constant(rho) && is_constant(cm)){struct{double x,y;}_const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann ((_const_alpha.x ? val(a.x,1,0,0)*val(fm.x,1,0,0)/_const_alpha.x : val(a.x,1,0,0)*_const_rho/(_const_cm + 1e-30)));}}}else if(!is_constant(alpha.x) && is_constant(a.x) && !is_constant(fm.x) && is_constant(rho) && is_constant(cm)){struct{double x,y;}_const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann ((val(alpha.x,1,0,0) ? _const_a.x*val(fm.x,1,0,0)/val(alpha.x,1,0,0) : _const_a.x*_const_rho/(_const_cm + 1e-30)));}}}else if(is_constant(alpha.x) && is_constant(a.x) && !is_constant(fm.x) && is_constant(rho) && is_constant(cm)){struct{double x,y;}_const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);struct{double x,y;}_const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann ((_const_alpha.x ? _const_a.x*val(fm.x,1,0,0)/_const_alpha.x : _const_a.x*_const_rho/(_const_cm + 1e-30)));}}}else if(!is_constant(alpha.x) && !is_constant(a.x) && is_constant(fm.x) && is_constant(rho) && is_constant(cm)){struct{double x,y;}_const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann ((val(alpha.x,1,0,0) ? val(a.x,1,0,0)*_const_fm.x/val(alpha.x,1,0,0) : val(a.x,1,0,0)*_const_rho/(_const_cm + 1e-30)));}}}else if(is_constant(alpha.x) && !is_constant(a.x) && is_constant(fm.x) && is_constant(rho) && is_constant(cm)){struct{double x,y;}_const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);struct{double x,y;}_const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann ((_const_alpha.x ? val(a.x,1,0,0)*_const_fm.x/_const_alpha.x : val(a.x,1,0,0)*_const_rho/(_const_cm + 1e-30)));}}}else if(!is_constant(alpha.x) && is_constant(a.x) && is_constant(fm.x) && is_constant(rho) && is_constant(cm)){struct{double x,y;}_const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);struct{double x,y;}_const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann ((val(alpha.x,1,0,0) ? _const_a.x*_const_fm.x/val(alpha.x,1,0,0) : _const_a.x*_const_rho/(_const_cm + 1e-30)));}}}else {struct{double x,y;}_const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);struct{double x,y;}_const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);struct{double x,y;}_const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann ((_const_alpha.x ? _const_a.x*_const_fm.x/_const_alpha.x : _const_a.x*_const_rho/(_const_cm + 1e-30)));}}}}static double _boundary0_homogeneous(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann_homogeneous ();}}
static double _boundary1(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;if(!is_constant(alpha.x) && !is_constant(a.x) && !is_constant(fm.x) && !is_constant(rho) && !is_constant(cm)){{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann (- (val(alpha.x,0,0,0) ? val(a.x,0,0,0)*val(fm.x,0,0,0)/val(alpha.x,0,0,0) : val(a.x,0,0,0)*val(rho,0,0,0)/(val(cm,0,0,0) + 1e-30)));}}}else if(is_constant(alpha.x) && !is_constant(a.x) && !is_constant(fm.x) && !is_constant(rho) && !is_constant(cm)){struct{double x,y;}_const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann (- (_const_alpha.x ? val(a.x,0,0,0)*val(fm.x,0,0,0)/_const_alpha.x : val(a.x,0,0,0)*val(rho,0,0,0)/(val(cm,0,0,0) + 1e-30)));}}}else if(!is_constant(alpha.x) && is_constant(a.x) && !is_constant(fm.x) && !is_constant(rho) && !is_constant(cm)){struct{double x,y;}_const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann (- (val(alpha.x,0,0,0) ? _const_a.x*val(fm.x,0,0,0)/val(alpha.x,0,0,0) : _const_a.x*val(rho,0,0,0)/(val(cm,0,0,0) + 1e-30)));}}}else if(is_constant(alpha.x) && is_constant(a.x) && !is_constant(fm.x) && !is_constant(rho) && !is_constant(cm)){struct{double x,y;}_const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);struct{double x,y;}_const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann (- (_const_alpha.x ? _const_a.x*val(fm.x,0,0,0)/_const_alpha.x : _const_a.x*val(rho,0,0,0)/(val(cm,0,0,0) + 1e-30)));}}}else if(!is_constant(alpha.x) && !is_constant(a.x) && is_constant(fm.x) && !is_constant(rho) && !is_constant(cm)){struct{double x,y;}_const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann (- (val(alpha.x,0,0,0) ? val(a.x,0,0,0)*_const_fm.x/val(alpha.x,0,0,0) : val(a.x,0,0,0)*val(rho,0,0,0)/(val(cm,0,0,0) + 1e-30)));}}}else if(is_constant(alpha.x) && !is_constant(a.x) && is_constant(fm.x) && !is_constant(rho) && !is_constant(cm)){struct{double x,y;}_const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);struct{double x,y;}_const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann (- (_const_alpha.x ? val(a.x,0,0,0)*_const_fm.x/_const_alpha.x : val(a.x,0,0,0)*val(rho,0,0,0)/(val(cm,0,0,0) + 1e-30)));}}}else if(!is_constant(alpha.x) && is_constant(a.x) && is_constant(fm.x) && !is_constant(rho) && !is_constant(cm)){struct{double x,y;}_const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);struct{double x,y;}_const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann (- (val(alpha.x,0,0,0) ? _const_a.x*_const_fm.x/val(alpha.x,0,0,0) : _const_a.x*val(rho,0,0,0)/(val(cm,0,0,0) + 1e-30)));}}}else if(is_constant(alpha.x) && is_constant(a.x) && is_constant(fm.x) && !is_constant(rho) && !is_constant(cm)){struct{double x,y;}_const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);struct{double x,y;}_const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);struct{double x,y;}_const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann (- (_const_alpha.x ? _const_a.x*_const_fm.x/_const_alpha.x : _const_a.x*val(rho,0,0,0)/(val(cm,0,0,0) + 1e-30)));}}}else if(!is_constant(alpha.x) && !is_constant(a.x) && !is_constant(fm.x) && is_constant(rho) && !is_constant(cm)){double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann (- (val(alpha.x,0,0,0) ? val(a.x,0,0,0)*val(fm.x,0,0,0)/val(alpha.x,0,0,0) : val(a.x,0,0,0)*_const_rho/(val(cm,0,0,0) + 1e-30)));}}}else if(is_constant(alpha.x) && !is_constant(a.x) && !is_constant(fm.x) && is_constant(rho) && !is_constant(cm)){struct{double x,y;}_const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann (- (_const_alpha.x ? val(a.x,0,0,0)*val(fm.x,0,0,0)/_const_alpha.x : val(a.x,0,0,0)*_const_rho/(val(cm,0,0,0) + 1e-30)));}}}else if(!is_constant(alpha.x) && is_constant(a.x) && !is_constant(fm.x) && is_constant(rho) && !is_constant(cm)){struct{double x,y;}_const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann (- (val(alpha.x,0,0,0) ? _const_a.x*val(fm.x,0,0,0)/val(alpha.x,0,0,0) : _const_a.x*_const_rho/(val(cm,0,0,0) + 1e-30)));}}}else if(is_constant(alpha.x) && is_constant(a.x) && !is_constant(fm.x) && is_constant(rho) && !is_constant(cm)){struct{double x,y;}_const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);struct{double x,y;}_const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann (- (_const_alpha.x ? _const_a.x*val(fm.x,0,0,0)/_const_alpha.x : _const_a.x*_const_rho/(val(cm,0,0,0) + 1e-30)));}}}else if(!is_constant(alpha.x) && !is_constant(a.x) && is_constant(fm.x) && is_constant(rho) && !is_constant(cm)){struct{double x,y;}_const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann (- (val(alpha.x,0,0,0) ? val(a.x,0,0,0)*_const_fm.x/val(alpha.x,0,0,0) : val(a.x,0,0,0)*_const_rho/(val(cm,0,0,0) + 1e-30)));}}}else if(is_constant(alpha.x) && !is_constant(a.x) && is_constant(fm.x) && is_constant(rho) && !is_constant(cm)){struct{double x,y;}_const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);struct{double x,y;}_const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann (- (_const_alpha.x ? val(a.x,0,0,0)*_const_fm.x/_const_alpha.x : val(a.x,0,0,0)*_const_rho/(val(cm,0,0,0) + 1e-30)));}}}else if(!is_constant(alpha.x) && is_constant(a.x) && is_constant(fm.x) && is_constant(rho) && !is_constant(cm)){struct{double x,y;}_const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);struct{double x,y;}_const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann (- (val(alpha.x,0,0,0) ? _const_a.x*_const_fm.x/val(alpha.x,0,0,0) : _const_a.x*_const_rho/(val(cm,0,0,0) + 1e-30)));}}}else if(is_constant(alpha.x) && is_constant(a.x) && is_constant(fm.x) && is_constant(rho) && !is_constant(cm)){struct{double x,y;}_const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);struct{double x,y;}_const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);struct{double x,y;}_const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann (- (_const_alpha.x ? _const_a.x*_const_fm.x/_const_alpha.x : _const_a.x*_const_rho/(val(cm,0,0,0) + 1e-30)));}}}else if(!is_constant(alpha.x) && !is_constant(a.x) && !is_constant(fm.x) && !is_constant(rho) && is_constant(cm)){double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann (- (val(alpha.x,0,0,0) ? val(a.x,0,0,0)*val(fm.x,0,0,0)/val(alpha.x,0,0,0) : val(a.x,0,0,0)*val(rho,0,0,0)/(_const_cm + 1e-30)));}}}else if(is_constant(alpha.x) && !is_constant(a.x) && !is_constant(fm.x) && !is_constant(rho) && is_constant(cm)){struct{double x,y;}_const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann (- (_const_alpha.x ? val(a.x,0,0,0)*val(fm.x,0,0,0)/_const_alpha.x : val(a.x,0,0,0)*val(rho,0,0,0)/(_const_cm + 1e-30)));}}}else if(!is_constant(alpha.x) && is_constant(a.x) && !is_constant(fm.x) && !is_constant(rho) && is_constant(cm)){struct{double x,y;}_const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann (- (val(alpha.x,0,0,0) ? _const_a.x*val(fm.x,0,0,0)/val(alpha.x,0,0,0) : _const_a.x*val(rho,0,0,0)/(_const_cm + 1e-30)));}}}else if(is_constant(alpha.x) && is_constant(a.x) && !is_constant(fm.x) && !is_constant(rho) && is_constant(cm)){struct{double x,y;}_const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);struct{double x,y;}_const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann (- (_const_alpha.x ? _const_a.x*val(fm.x,0,0,0)/_const_alpha.x : _const_a.x*val(rho,0,0,0)/(_const_cm + 1e-30)));}}}else if(!is_constant(alpha.x) && !is_constant(a.x) && is_constant(fm.x) && !is_constant(rho) && is_constant(cm)){struct{double x,y;}_const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann (- (val(alpha.x,0,0,0) ? val(a.x,0,0,0)*_const_fm.x/val(alpha.x,0,0,0) : val(a.x,0,0,0)*val(rho,0,0,0)/(_const_cm + 1e-30)));}}}else if(is_constant(alpha.x) && !is_constant(a.x) && is_constant(fm.x) && !is_constant(rho) && is_constant(cm)){struct{double x,y;}_const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);struct{double x,y;}_const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann (- (_const_alpha.x ? val(a.x,0,0,0)*_const_fm.x/_const_alpha.x : val(a.x,0,0,0)*val(rho,0,0,0)/(_const_cm + 1e-30)));}}}else if(!is_constant(alpha.x) && is_constant(a.x) && is_constant(fm.x) && !is_constant(rho) && is_constant(cm)){struct{double x,y;}_const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);struct{double x,y;}_const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann (- (val(alpha.x,0,0,0) ? _const_a.x*_const_fm.x/val(alpha.x,0,0,0) : _const_a.x*val(rho,0,0,0)/(_const_cm + 1e-30)));}}}else if(is_constant(alpha.x) && is_constant(a.x) && is_constant(fm.x) && !is_constant(rho) && is_constant(cm)){struct{double x,y;}_const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);struct{double x,y;}_const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);struct{double x,y;}_const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann (- (_const_alpha.x ? _const_a.x*_const_fm.x/_const_alpha.x : _const_a.x*val(rho,0,0,0)/(_const_cm + 1e-30)));}}}else if(!is_constant(alpha.x) && !is_constant(a.x) && !is_constant(fm.x) && is_constant(rho) && is_constant(cm)){double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann (- (val(alpha.x,0,0,0) ? val(a.x,0,0,0)*val(fm.x,0,0,0)/val(alpha.x,0,0,0) : val(a.x,0,0,0)*_const_rho/(_const_cm + 1e-30)));}}}else if(is_constant(alpha.x) && !is_constant(a.x) && !is_constant(fm.x) && is_constant(rho) && is_constant(cm)){struct{double x,y;}_const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann (- (_const_alpha.x ? val(a.x,0,0,0)*val(fm.x,0,0,0)/_const_alpha.x : val(a.x,0,0,0)*_const_rho/(_const_cm + 1e-30)));}}}else if(!is_constant(alpha.x) && is_constant(a.x) && !is_constant(fm.x) && is_constant(rho) && is_constant(cm)){struct{double x,y;}_const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann (- (val(alpha.x,0,0,0) ? _const_a.x*val(fm.x,0,0,0)/val(alpha.x,0,0,0) : _const_a.x*_const_rho/(_const_cm + 1e-30)));}}}else if(is_constant(alpha.x) && is_constant(a.x) && !is_constant(fm.x) && is_constant(rho) && is_constant(cm)){struct{double x,y;}_const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);struct{double x,y;}_const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann (- (_const_alpha.x ? _const_a.x*val(fm.x,0,0,0)/_const_alpha.x : _const_a.x*_const_rho/(_const_cm + 1e-30)));}}}else if(!is_constant(alpha.x) && !is_constant(a.x) && is_constant(fm.x) && is_constant(rho) && is_constant(cm)){struct{double x,y;}_const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann (- (val(alpha.x,0,0,0) ? val(a.x,0,0,0)*_const_fm.x/val(alpha.x,0,0,0) : val(a.x,0,0,0)*_const_rho/(_const_cm + 1e-30)));}}}else if(is_constant(alpha.x) && !is_constant(a.x) && is_constant(fm.x) && is_constant(rho) && is_constant(cm)){struct{double x,y;}_const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);struct{double x,y;}_const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann (- (_const_alpha.x ? val(a.x,0,0,0)*_const_fm.x/_const_alpha.x : val(a.x,0,0,0)*_const_rho/(_const_cm + 1e-30)));}}}else if(!is_constant(alpha.x) && is_constant(a.x) && is_constant(fm.x) && is_constant(rho) && is_constant(cm)){struct{double x,y;}_const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);struct{double x,y;}_const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann (- (val(alpha.x,0,0,0) ? _const_a.x*_const_fm.x/val(alpha.x,0,0,0) : _const_a.x*_const_rho/(_const_cm + 1e-30)));}}}else {struct{double x,y;}_const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);struct{double x,y;}_const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);struct{double x,y;}_const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann (- (_const_alpha.x ? _const_a.x*_const_fm.x/_const_alpha.x : _const_a.x*_const_rho/(_const_cm + 1e-30)));}}}}static double _boundary1_homogeneous(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann_homogeneous ();}}








static double _boundary2(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;if(!is_constant(alpha.x) && !is_constant(a.x) && !is_constant(fm.x) && !is_constant(rho) && !is_constant(cm)){{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann ((val(alpha.y,0,1,0) ? val(a.y,0,1,0)*val(fm.y,0,1,0)/val(alpha.y,0,1,0) : val(a.y,0,1,0)*val(rho,0,0,0)/(val(cm,0,0,0) + 1e-30)));}}}else if(is_constant(alpha.x) && !is_constant(a.x) && !is_constant(fm.x) && !is_constant(rho) && !is_constant(cm)){struct{double x,y;}_const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann ((_const_alpha.y ? val(a.y,0,1,0)*val(fm.y,0,1,0)/_const_alpha.y : val(a.y,0,1,0)*val(rho,0,0,0)/(val(cm,0,0,0) + 1e-30)));}}}else if(!is_constant(alpha.x) && is_constant(a.x) && !is_constant(fm.x) && !is_constant(rho) && !is_constant(cm)){struct{double x,y;}_const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann ((val(alpha.y,0,1,0) ? _const_a.y*val(fm.y,0,1,0)/val(alpha.y,0,1,0) : _const_a.y*val(rho,0,0,0)/(val(cm,0,0,0) + 1e-30)));}}}else if(is_constant(alpha.x) && is_constant(a.x) && !is_constant(fm.x) && !is_constant(rho) && !is_constant(cm)){struct{double x,y;}_const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);struct{double x,y;}_const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann ((_const_alpha.y ? _const_a.y*val(fm.y,0,1,0)/_const_alpha.y : _const_a.y*val(rho,0,0,0)/(val(cm,0,0,0) + 1e-30)));}}}else if(!is_constant(alpha.x) && !is_constant(a.x) && is_constant(fm.x) && !is_constant(rho) && !is_constant(cm)){struct{double x,y;}_const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann ((val(alpha.y,0,1,0) ? val(a.y,0,1,0)*_const_fm.y/val(alpha.y,0,1,0) : val(a.y,0,1,0)*val(rho,0,0,0)/(val(cm,0,0,0) + 1e-30)));}}}else if(is_constant(alpha.x) && !is_constant(a.x) && is_constant(fm.x) && !is_constant(rho) && !is_constant(cm)){struct{double x,y;}_const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);struct{double x,y;}_const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann ((_const_alpha.y ? val(a.y,0,1,0)*_const_fm.y/_const_alpha.y : val(a.y,0,1,0)*val(rho,0,0,0)/(val(cm,0,0,0) + 1e-30)));}}}else if(!is_constant(alpha.x) && is_constant(a.x) && is_constant(fm.x) && !is_constant(rho) && !is_constant(cm)){struct{double x,y;}_const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);struct{double x,y;}_const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann ((val(alpha.y,0,1,0) ? _const_a.y*_const_fm.y/val(alpha.y,0,1,0) : _const_a.y*val(rho,0,0,0)/(val(cm,0,0,0) + 1e-30)));}}}else if(is_constant(alpha.x) && is_constant(a.x) && is_constant(fm.x) && !is_constant(rho) && !is_constant(cm)){struct{double x,y;}_const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);struct{double x,y;}_const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);struct{double x,y;}_const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann ((_const_alpha.y ? _const_a.y*_const_fm.y/_const_alpha.y : _const_a.y*val(rho,0,0,0)/(val(cm,0,0,0) + 1e-30)));}}}else if(!is_constant(alpha.x) && !is_constant(a.x) && !is_constant(fm.x) && is_constant(rho) && !is_constant(cm)){double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann ((val(alpha.y,0,1,0) ? val(a.y,0,1,0)*val(fm.y,0,1,0)/val(alpha.y,0,1,0) : val(a.y,0,1,0)*_const_rho/(val(cm,0,0,0) + 1e-30)));}}}else if(is_constant(alpha.x) && !is_constant(a.x) && !is_constant(fm.x) && is_constant(rho) && !is_constant(cm)){struct{double x,y;}_const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann ((_const_alpha.y ? val(a.y,0,1,0)*val(fm.y,0,1,0)/_const_alpha.y : val(a.y,0,1,0)*_const_rho/(val(cm,0,0,0) + 1e-30)));}}}else if(!is_constant(alpha.x) && is_constant(a.x) && !is_constant(fm.x) && is_constant(rho) && !is_constant(cm)){struct{double x,y;}_const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann ((val(alpha.y,0,1,0) ? _const_a.y*val(fm.y,0,1,0)/val(alpha.y,0,1,0) : _const_a.y*_const_rho/(val(cm,0,0,0) + 1e-30)));}}}else if(is_constant(alpha.x) && is_constant(a.x) && !is_constant(fm.x) && is_constant(rho) && !is_constant(cm)){struct{double x,y;}_const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);struct{double x,y;}_const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann ((_const_alpha.y ? _const_a.y*val(fm.y,0,1,0)/_const_alpha.y : _const_a.y*_const_rho/(val(cm,0,0,0) + 1e-30)));}}}else if(!is_constant(alpha.x) && !is_constant(a.x) && is_constant(fm.x) && is_constant(rho) && !is_constant(cm)){struct{double x,y;}_const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann ((val(alpha.y,0,1,0) ? val(a.y,0,1,0)*_const_fm.y/val(alpha.y,0,1,0) : val(a.y,0,1,0)*_const_rho/(val(cm,0,0,0) + 1e-30)));}}}else if(is_constant(alpha.x) && !is_constant(a.x) && is_constant(fm.x) && is_constant(rho) && !is_constant(cm)){struct{double x,y;}_const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);struct{double x,y;}_const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann ((_const_alpha.y ? val(a.y,0,1,0)*_const_fm.y/_const_alpha.y : val(a.y,0,1,0)*_const_rho/(val(cm,0,0,0) + 1e-30)));}}}else if(!is_constant(alpha.x) && is_constant(a.x) && is_constant(fm.x) && is_constant(rho) && !is_constant(cm)){struct{double x,y;}_const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);struct{double x,y;}_const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann ((val(alpha.y,0,1,0) ? _const_a.y*_const_fm.y/val(alpha.y,0,1,0) : _const_a.y*_const_rho/(val(cm,0,0,0) + 1e-30)));}}}else if(is_constant(alpha.x) && is_constant(a.x) && is_constant(fm.x) && is_constant(rho) && !is_constant(cm)){struct{double x,y;}_const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);struct{double x,y;}_const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);struct{double x,y;}_const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann ((_const_alpha.y ? _const_a.y*_const_fm.y/_const_alpha.y : _const_a.y*_const_rho/(val(cm,0,0,0) + 1e-30)));}}}else if(!is_constant(alpha.x) && !is_constant(a.x) && !is_constant(fm.x) && !is_constant(rho) && is_constant(cm)){double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann ((val(alpha.y,0,1,0) ? val(a.y,0,1,0)*val(fm.y,0,1,0)/val(alpha.y,0,1,0) : val(a.y,0,1,0)*val(rho,0,0,0)/(_const_cm + 1e-30)));}}}else if(is_constant(alpha.x) && !is_constant(a.x) && !is_constant(fm.x) && !is_constant(rho) && is_constant(cm)){struct{double x,y;}_const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann ((_const_alpha.y ? val(a.y,0,1,0)*val(fm.y,0,1,0)/_const_alpha.y : val(a.y,0,1,0)*val(rho,0,0,0)/(_const_cm + 1e-30)));}}}else if(!is_constant(alpha.x) && is_constant(a.x) && !is_constant(fm.x) && !is_constant(rho) && is_constant(cm)){struct{double x,y;}_const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann ((val(alpha.y,0,1,0) ? _const_a.y*val(fm.y,0,1,0)/val(alpha.y,0,1,0) : _const_a.y*val(rho,0,0,0)/(_const_cm + 1e-30)));}}}else if(is_constant(alpha.x) && is_constant(a.x) && !is_constant(fm.x) && !is_constant(rho) && is_constant(cm)){struct{double x,y;}_const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);struct{double x,y;}_const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann ((_const_alpha.y ? _const_a.y*val(fm.y,0,1,0)/_const_alpha.y : _const_a.y*val(rho,0,0,0)/(_const_cm + 1e-30)));}}}else if(!is_constant(alpha.x) && !is_constant(a.x) && is_constant(fm.x) && !is_constant(rho) && is_constant(cm)){struct{double x,y;}_const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann ((val(alpha.y,0,1,0) ? val(a.y,0,1,0)*_const_fm.y/val(alpha.y,0,1,0) : val(a.y,0,1,0)*val(rho,0,0,0)/(_const_cm + 1e-30)));}}}else if(is_constant(alpha.x) && !is_constant(a.x) && is_constant(fm.x) && !is_constant(rho) && is_constant(cm)){struct{double x,y;}_const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);struct{double x,y;}_const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann ((_const_alpha.y ? val(a.y,0,1,0)*_const_fm.y/_const_alpha.y : val(a.y,0,1,0)*val(rho,0,0,0)/(_const_cm + 1e-30)));}}}else if(!is_constant(alpha.x) && is_constant(a.x) && is_constant(fm.x) && !is_constant(rho) && is_constant(cm)){struct{double x,y;}_const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);struct{double x,y;}_const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann ((val(alpha.y,0,1,0) ? _const_a.y*_const_fm.y/val(alpha.y,0,1,0) : _const_a.y*val(rho,0,0,0)/(_const_cm + 1e-30)));}}}else if(is_constant(alpha.x) && is_constant(a.x) && is_constant(fm.x) && !is_constant(rho) && is_constant(cm)){struct{double x,y;}_const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);struct{double x,y;}_const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);struct{double x,y;}_const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann ((_const_alpha.y ? _const_a.y*_const_fm.y/_const_alpha.y : _const_a.y*val(rho,0,0,0)/(_const_cm + 1e-30)));}}}else if(!is_constant(alpha.x) && !is_constant(a.x) && !is_constant(fm.x) && is_constant(rho) && is_constant(cm)){double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann ((val(alpha.y,0,1,0) ? val(a.y,0,1,0)*val(fm.y,0,1,0)/val(alpha.y,0,1,0) : val(a.y,0,1,0)*_const_rho/(_const_cm + 1e-30)));}}}else if(is_constant(alpha.x) && !is_constant(a.x) && !is_constant(fm.x) && is_constant(rho) && is_constant(cm)){struct{double x,y;}_const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann ((_const_alpha.y ? val(a.y,0,1,0)*val(fm.y,0,1,0)/_const_alpha.y : val(a.y,0,1,0)*_const_rho/(_const_cm + 1e-30)));}}}else if(!is_constant(alpha.x) && is_constant(a.x) && !is_constant(fm.x) && is_constant(rho) && is_constant(cm)){struct{double x,y;}_const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann ((val(alpha.y,0,1,0) ? _const_a.y*val(fm.y,0,1,0)/val(alpha.y,0,1,0) : _const_a.y*_const_rho/(_const_cm + 1e-30)));}}}else if(is_constant(alpha.x) && is_constant(a.x) && !is_constant(fm.x) && is_constant(rho) && is_constant(cm)){struct{double x,y;}_const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);struct{double x,y;}_const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann ((_const_alpha.y ? _const_a.y*val(fm.y,0,1,0)/_const_alpha.y : _const_a.y*_const_rho/(_const_cm + 1e-30)));}}}else if(!is_constant(alpha.x) && !is_constant(a.x) && is_constant(fm.x) && is_constant(rho) && is_constant(cm)){struct{double x,y;}_const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann ((val(alpha.y,0,1,0) ? val(a.y,0,1,0)*_const_fm.y/val(alpha.y,0,1,0) : val(a.y,0,1,0)*_const_rho/(_const_cm + 1e-30)));}}}else if(is_constant(alpha.x) && !is_constant(a.x) && is_constant(fm.x) && is_constant(rho) && is_constant(cm)){struct{double x,y;}_const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);struct{double x,y;}_const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann ((_const_alpha.y ? val(a.y,0,1,0)*_const_fm.y/_const_alpha.y : val(a.y,0,1,0)*_const_rho/(_const_cm + 1e-30)));}}}else if(!is_constant(alpha.x) && is_constant(a.x) && is_constant(fm.x) && is_constant(rho) && is_constant(cm)){struct{double x,y;}_const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);struct{double x,y;}_const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann ((val(alpha.y,0,1,0) ? _const_a.y*_const_fm.y/val(alpha.y,0,1,0) : _const_a.y*_const_rho/(_const_cm + 1e-30)));}}}else {struct{double x,y;}_const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);struct{double x,y;}_const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);struct{double x,y;}_const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann ((_const_alpha.y ? _const_a.y*_const_fm.y/_const_alpha.y : _const_a.y*_const_rho/(_const_cm + 1e-30)));}}}}static double _boundary2_homogeneous(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann_homogeneous ();}}
static double _boundary3(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;if(!is_constant(alpha.x) && !is_constant(a.x) && !is_constant(fm.x) && !is_constant(rho) && !is_constant(cm)){{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann (- (val(alpha.y,0,0,0) ? val(a.y,0,0,0)*val(fm.y,0,0,0)/val(alpha.y,0,0,0) : val(a.y,0,0,0)*val(rho,0,0,0)/(val(cm,0,0,0) + 1e-30)));}}}else if(is_constant(alpha.x) && !is_constant(a.x) && !is_constant(fm.x) && !is_constant(rho) && !is_constant(cm)){struct{double x,y;}_const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann (- (_const_alpha.y ? val(a.y,0,0,0)*val(fm.y,0,0,0)/_const_alpha.y : val(a.y,0,0,0)*val(rho,0,0,0)/(val(cm,0,0,0) + 1e-30)));}}}else if(!is_constant(alpha.x) && is_constant(a.x) && !is_constant(fm.x) && !is_constant(rho) && !is_constant(cm)){struct{double x,y;}_const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann (- (val(alpha.y,0,0,0) ? _const_a.y*val(fm.y,0,0,0)/val(alpha.y,0,0,0) : _const_a.y*val(rho,0,0,0)/(val(cm,0,0,0) + 1e-30)));}}}else if(is_constant(alpha.x) && is_constant(a.x) && !is_constant(fm.x) && !is_constant(rho) && !is_constant(cm)){struct{double x,y;}_const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);struct{double x,y;}_const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann (- (_const_alpha.y ? _const_a.y*val(fm.y,0,0,0)/_const_alpha.y : _const_a.y*val(rho,0,0,0)/(val(cm,0,0,0) + 1e-30)));}}}else if(!is_constant(alpha.x) && !is_constant(a.x) && is_constant(fm.x) && !is_constant(rho) && !is_constant(cm)){struct{double x,y;}_const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann (- (val(alpha.y,0,0,0) ? val(a.y,0,0,0)*_const_fm.y/val(alpha.y,0,0,0) : val(a.y,0,0,0)*val(rho,0,0,0)/(val(cm,0,0,0) + 1e-30)));}}}else if(is_constant(alpha.x) && !is_constant(a.x) && is_constant(fm.x) && !is_constant(rho) && !is_constant(cm)){struct{double x,y;}_const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);struct{double x,y;}_const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann (- (_const_alpha.y ? val(a.y,0,0,0)*_const_fm.y/_const_alpha.y : val(a.y,0,0,0)*val(rho,0,0,0)/(val(cm,0,0,0) + 1e-30)));}}}else if(!is_constant(alpha.x) && is_constant(a.x) && is_constant(fm.x) && !is_constant(rho) && !is_constant(cm)){struct{double x,y;}_const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);struct{double x,y;}_const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann (- (val(alpha.y,0,0,0) ? _const_a.y*_const_fm.y/val(alpha.y,0,0,0) : _const_a.y*val(rho,0,0,0)/(val(cm,0,0,0) + 1e-30)));}}}else if(is_constant(alpha.x) && is_constant(a.x) && is_constant(fm.x) && !is_constant(rho) && !is_constant(cm)){struct{double x,y;}_const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);struct{double x,y;}_const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);struct{double x,y;}_const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann (- (_const_alpha.y ? _const_a.y*_const_fm.y/_const_alpha.y : _const_a.y*val(rho,0,0,0)/(val(cm,0,0,0) + 1e-30)));}}}else if(!is_constant(alpha.x) && !is_constant(a.x) && !is_constant(fm.x) && is_constant(rho) && !is_constant(cm)){double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann (- (val(alpha.y,0,0,0) ? val(a.y,0,0,0)*val(fm.y,0,0,0)/val(alpha.y,0,0,0) : val(a.y,0,0,0)*_const_rho/(val(cm,0,0,0) + 1e-30)));}}}else if(is_constant(alpha.x) && !is_constant(a.x) && !is_constant(fm.x) && is_constant(rho) && !is_constant(cm)){struct{double x,y;}_const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann (- (_const_alpha.y ? val(a.y,0,0,0)*val(fm.y,0,0,0)/_const_alpha.y : val(a.y,0,0,0)*_const_rho/(val(cm,0,0,0) + 1e-30)));}}}else if(!is_constant(alpha.x) && is_constant(a.x) && !is_constant(fm.x) && is_constant(rho) && !is_constant(cm)){struct{double x,y;}_const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann (- (val(alpha.y,0,0,0) ? _const_a.y*val(fm.y,0,0,0)/val(alpha.y,0,0,0) : _const_a.y*_const_rho/(val(cm,0,0,0) + 1e-30)));}}}else if(is_constant(alpha.x) && is_constant(a.x) && !is_constant(fm.x) && is_constant(rho) && !is_constant(cm)){struct{double x,y;}_const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);struct{double x,y;}_const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann (- (_const_alpha.y ? _const_a.y*val(fm.y,0,0,0)/_const_alpha.y : _const_a.y*_const_rho/(val(cm,0,0,0) + 1e-30)));}}}else if(!is_constant(alpha.x) && !is_constant(a.x) && is_constant(fm.x) && is_constant(rho) && !is_constant(cm)){struct{double x,y;}_const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann (- (val(alpha.y,0,0,0) ? val(a.y,0,0,0)*_const_fm.y/val(alpha.y,0,0,0) : val(a.y,0,0,0)*_const_rho/(val(cm,0,0,0) + 1e-30)));}}}else if(is_constant(alpha.x) && !is_constant(a.x) && is_constant(fm.x) && is_constant(rho) && !is_constant(cm)){struct{double x,y;}_const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);struct{double x,y;}_const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann (- (_const_alpha.y ? val(a.y,0,0,0)*_const_fm.y/_const_alpha.y : val(a.y,0,0,0)*_const_rho/(val(cm,0,0,0) + 1e-30)));}}}else if(!is_constant(alpha.x) && is_constant(a.x) && is_constant(fm.x) && is_constant(rho) && !is_constant(cm)){struct{double x,y;}_const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);struct{double x,y;}_const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann (- (val(alpha.y,0,0,0) ? _const_a.y*_const_fm.y/val(alpha.y,0,0,0) : _const_a.y*_const_rho/(val(cm,0,0,0) + 1e-30)));}}}else if(is_constant(alpha.x) && is_constant(a.x) && is_constant(fm.x) && is_constant(rho) && !is_constant(cm)){struct{double x,y;}_const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);struct{double x,y;}_const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);struct{double x,y;}_const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann (- (_const_alpha.y ? _const_a.y*_const_fm.y/_const_alpha.y : _const_a.y*_const_rho/(val(cm,0,0,0) + 1e-30)));}}}else if(!is_constant(alpha.x) && !is_constant(a.x) && !is_constant(fm.x) && !is_constant(rho) && is_constant(cm)){double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann (- (val(alpha.y,0,0,0) ? val(a.y,0,0,0)*val(fm.y,0,0,0)/val(alpha.y,0,0,0) : val(a.y,0,0,0)*val(rho,0,0,0)/(_const_cm + 1e-30)));}}}else if(is_constant(alpha.x) && !is_constant(a.x) && !is_constant(fm.x) && !is_constant(rho) && is_constant(cm)){struct{double x,y;}_const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann (- (_const_alpha.y ? val(a.y,0,0,0)*val(fm.y,0,0,0)/_const_alpha.y : val(a.y,0,0,0)*val(rho,0,0,0)/(_const_cm + 1e-30)));}}}else if(!is_constant(alpha.x) && is_constant(a.x) && !is_constant(fm.x) && !is_constant(rho) && is_constant(cm)){struct{double x,y;}_const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann (- (val(alpha.y,0,0,0) ? _const_a.y*val(fm.y,0,0,0)/val(alpha.y,0,0,0) : _const_a.y*val(rho,0,0,0)/(_const_cm + 1e-30)));}}}else if(is_constant(alpha.x) && is_constant(a.x) && !is_constant(fm.x) && !is_constant(rho) && is_constant(cm)){struct{double x,y;}_const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);struct{double x,y;}_const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann (- (_const_alpha.y ? _const_a.y*val(fm.y,0,0,0)/_const_alpha.y : _const_a.y*val(rho,0,0,0)/(_const_cm + 1e-30)));}}}else if(!is_constant(alpha.x) && !is_constant(a.x) && is_constant(fm.x) && !is_constant(rho) && is_constant(cm)){struct{double x,y;}_const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann (- (val(alpha.y,0,0,0) ? val(a.y,0,0,0)*_const_fm.y/val(alpha.y,0,0,0) : val(a.y,0,0,0)*val(rho,0,0,0)/(_const_cm + 1e-30)));}}}else if(is_constant(alpha.x) && !is_constant(a.x) && is_constant(fm.x) && !is_constant(rho) && is_constant(cm)){struct{double x,y;}_const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);struct{double x,y;}_const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann (- (_const_alpha.y ? val(a.y,0,0,0)*_const_fm.y/_const_alpha.y : val(a.y,0,0,0)*val(rho,0,0,0)/(_const_cm + 1e-30)));}}}else if(!is_constant(alpha.x) && is_constant(a.x) && is_constant(fm.x) && !is_constant(rho) && is_constant(cm)){struct{double x,y;}_const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);struct{double x,y;}_const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann (- (val(alpha.y,0,0,0) ? _const_a.y*_const_fm.y/val(alpha.y,0,0,0) : _const_a.y*val(rho,0,0,0)/(_const_cm + 1e-30)));}}}else if(is_constant(alpha.x) && is_constant(a.x) && is_constant(fm.x) && !is_constant(rho) && is_constant(cm)){struct{double x,y;}_const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);struct{double x,y;}_const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);struct{double x,y;}_const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann (- (_const_alpha.y ? _const_a.y*_const_fm.y/_const_alpha.y : _const_a.y*val(rho,0,0,0)/(_const_cm + 1e-30)));}}}else if(!is_constant(alpha.x) && !is_constant(a.x) && !is_constant(fm.x) && is_constant(rho) && is_constant(cm)){double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann (- (val(alpha.y,0,0,0) ? val(a.y,0,0,0)*val(fm.y,0,0,0)/val(alpha.y,0,0,0) : val(a.y,0,0,0)*_const_rho/(_const_cm + 1e-30)));}}}else if(is_constant(alpha.x) && !is_constant(a.x) && !is_constant(fm.x) && is_constant(rho) && is_constant(cm)){struct{double x,y;}_const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann (- (_const_alpha.y ? val(a.y,0,0,0)*val(fm.y,0,0,0)/_const_alpha.y : val(a.y,0,0,0)*_const_rho/(_const_cm + 1e-30)));}}}else if(!is_constant(alpha.x) && is_constant(a.x) && !is_constant(fm.x) && is_constant(rho) && is_constant(cm)){struct{double x,y;}_const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann (- (val(alpha.y,0,0,0) ? _const_a.y*val(fm.y,0,0,0)/val(alpha.y,0,0,0) : _const_a.y*_const_rho/(_const_cm + 1e-30)));}}}else if(is_constant(alpha.x) && is_constant(a.x) && !is_constant(fm.x) && is_constant(rho) && is_constant(cm)){struct{double x,y;}_const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);struct{double x,y;}_const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann (- (_const_alpha.y ? _const_a.y*val(fm.y,0,0,0)/_const_alpha.y : _const_a.y*_const_rho/(_const_cm + 1e-30)));}}}else if(!is_constant(alpha.x) && !is_constant(a.x) && is_constant(fm.x) && is_constant(rho) && is_constant(cm)){struct{double x,y;}_const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann (- (val(alpha.y,0,0,0) ? val(a.y,0,0,0)*_const_fm.y/val(alpha.y,0,0,0) : val(a.y,0,0,0)*_const_rho/(_const_cm + 1e-30)));}}}else if(is_constant(alpha.x) && !is_constant(a.x) && is_constant(fm.x) && is_constant(rho) && is_constant(cm)){struct{double x,y;}_const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);struct{double x,y;}_const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann (- (_const_alpha.y ? val(a.y,0,0,0)*_const_fm.y/_const_alpha.y : val(a.y,0,0,0)*_const_rho/(_const_cm + 1e-30)));}}}else if(!is_constant(alpha.x) && is_constant(a.x) && is_constant(fm.x) && is_constant(rho) && is_constant(cm)){struct{double x,y;}_const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);struct{double x,y;}_const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann (- (val(alpha.y,0,0,0) ? _const_a.y*_const_fm.y/val(alpha.y,0,0,0) : _const_a.y*_const_rho/(_const_cm + 1e-30)));}}}else {struct{double x,y;}_const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);struct{double x,y;}_const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);struct{double x,y;}_const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann (- (_const_alpha.y ? _const_a.y*_const_fm.y/_const_alpha.y : _const_a.y*_const_rho/(_const_cm + 1e-30)));}}}}static double _boundary3_homogeneous(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann_homogeneous ();}}
#line 116 "/u/basilisk/src/navier-stokes/centered.h"
void pressure_embed_gradient (Point point, scalar p, coord * g)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;

#line 117
if(!is_constant(rho) && !is_constant(cm) && !is_constant(a.x)){{
  
    g->x = val(rho,0,0,0)/(val(cm,0,0,0) + 1e-30)*(val(a.x,0,0,0) + val(a.x,1,0,0))/2.;
    
#line 119
g->y = val(rho,0,0,0)/(val(cm,0,0,0) + 1e-30)*(val(a.y,0,0,0) + val(a.y,0,1,0))/2.;
}}else if(is_constant(rho) && !is_constant(cm) && !is_constant(a.x)){double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);

#line 117
{
  
    g->x = _const_rho/(val(cm,0,0,0) + 1e-30)*(val(a.x,0,0,0) + val(a.x,1,0,0))/2.;
    
#line 119
g->y = _const_rho/(val(cm,0,0,0) + 1e-30)*(val(a.y,0,0,0) + val(a.y,0,1,0))/2.;
}}else if(!is_constant(rho) && is_constant(cm) && !is_constant(a.x)){double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);

#line 117
{
  
    g->x = val(rho,0,0,0)/(_const_cm + 1e-30)*(val(a.x,0,0,0) + val(a.x,1,0,0))/2.;
    
#line 119
g->y = val(rho,0,0,0)/(_const_cm + 1e-30)*(val(a.y,0,0,0) + val(a.y,0,1,0))/2.;
}}else if(is_constant(rho) && is_constant(cm) && !is_constant(a.x)){double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);

#line 117
{
  
    g->x = _const_rho/(_const_cm + 1e-30)*(val(a.x,0,0,0) + val(a.x,1,0,0))/2.;
    
#line 119
g->y = _const_rho/(_const_cm + 1e-30)*(val(a.y,0,0,0) + val(a.y,0,1,0))/2.;
}}else if(!is_constant(rho) && !is_constant(cm) && is_constant(a.x)){struct{double x,y;}_const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);

#line 117
{
  
    g->x = val(rho,0,0,0)/(val(cm,0,0,0) + 1e-30)*(_const_a.x + _const_a.x)/2.;
    
#line 119
g->y = val(rho,0,0,0)/(val(cm,0,0,0) + 1e-30)*(_const_a.y + _const_a.y)/2.;
}}else if(is_constant(rho) && !is_constant(cm) && is_constant(a.x)){double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);struct{double x,y;}_const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);

#line 117
{
  
    g->x = _const_rho/(val(cm,0,0,0) + 1e-30)*(_const_a.x + _const_a.x)/2.;
    
#line 119
g->y = _const_rho/(val(cm,0,0,0) + 1e-30)*(_const_a.y + _const_a.y)/2.;
}}else if(!is_constant(rho) && is_constant(cm) && is_constant(a.x)){double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);struct{double x,y;}_const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);

#line 117
{
  
    g->x = val(rho,0,0,0)/(_const_cm + 1e-30)*(_const_a.x + _const_a.x)/2.;
    
#line 119
g->y = val(rho,0,0,0)/(_const_cm + 1e-30)*(_const_a.y + _const_a.y)/2.;
}}else {double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);struct{double x,y;}_const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);

#line 117
{
  
    g->x = _const_rho/(_const_cm + 1e-30)*(_const_a.x + _const_a.x)/2.;
    
#line 119
g->y = _const_rho/(_const_cm + 1e-30)*(_const_a.y + _const_a.y)/2.;
}}

#line 120
}





static int defaults_1_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i = 0);*ip=i;*tp=t;return ret;}      static int defaults_1(const int i,const double t,Event *_ev){tracing("defaults_1","/u/basilisk/src/navier-stokes/centered.h",0);
{

  CFL = 0.8;




  _attribute[p.i].nodump = _attribute[pf.i].nodump = true;




  if (alpha.x.i == unityf.x.i) {
    alpha = fm;
    rho = cm;
  }
  else if (!is_constant(alpha.x)) {
    vector alphav = alpha;
    foreach_face_stencil(){_stencil_is_face_x(){
      {_stencil_val_a(alphav.x,0,0,0); _stencil_val(fm.x,0,0,0); }}end__stencil_is_face_x()
#line 145
_stencil_is_face_y(){
      {_stencil_val_a(alphav.y,0,0,0); _stencil_val(fm.y,0,0,0); }}end__stencil_is_face_y()}end_foreach_face_stencil();
    
#line 145
if(!is_constant(fm.x)){{foreach_face_generic(){is_face_x(){
      val(alphav.x,0,0,0) = val(fm.x,0,0,0);}end_is_face_x()
#line 145
is_face_y(){
      val(alphav.y,0,0,0) = val(fm.y,0,0,0);}end_is_face_y()}end_foreach_face_generic();}}else {struct{double x,y;}_const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);
    {
#line 145
foreach_face_generic(){is_face_x(){
      val(alphav.x,0,0,0) = _const_fm.x;}end_is_face_x()
#line 145
is_face_y(){
      val(alphav.y,0,0,0) = _const_fm.y;}end_is_face_y()}end_foreach_face_generic();}}
  }






  _attribute[uf.x.i].refine = refine_face_solenoidal;






  _attribute[uf.x.i].refine = refine_face;
  
    _attribute[uf.x.i].prolongation = refine_embed_face_x;
    
#line 163
_attribute[uf.y.i].prolongation = refine_embed_face_y;
  {scalar*_i=(scalar*) ((scalar[]){p, pf, u.x,u.y, g.x,g.y,{-1}});if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {
    _attribute[s.i].restriction = restriction_embed_linear;
    _attribute[s.i].refine = _attribute[s.i].prolongation = refine_embed_linear;
    _attribute[s.i].depends = list_add (_attribute[s.i].depends, cs);
  }}}
  {scalar*_i=(scalar*) ((scalar[]){p, pf,{-1}});if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    _attribute[s.i].embed_gradient = pressure_embed_gradient;}}


}{end_tracing("defaults_1","/u/basilisk/src/navier-stokes/centered.h",0);return 0;}end_tracing("defaults_1","/u/basilisk/src/navier-stokes/centered.h",0);}





static int default_display_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i = 0);*ip=i;*tp=t;return ret;}      static int default_display(const int i,const double t,Event *_ev){tracing("default_display","/u/basilisk/src/navier-stokes/centered.h",0);
  display ((struct _display){"squares (color = 'u.x', spread = -1);"});{end_tracing("default_display","/u/basilisk/src/navier-stokes/centered.h",0);return 0;}end_tracing("default_display","/u/basilisk/src/navier-stokes/centered.h",0);}





double dtmax;

static int init_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i = 0);*ip=i;*tp=t;return ret;}      static int init(const int i,const double t,Event *_ev){tracing("init","/u/basilisk/src/navier-stokes/centered.h",0);
{
  trash (((vector[]){uf,{{-1},{-1}}}));
  foreach_face_stencil(){_stencil_is_face_x(){
    {_stencil_val_a(uf.x,0,0,0); _stencil_val(fm.x,0,0,0); _stencil_val(fs.x,0,0,0); _stencil_val(fs.x,0,0,0); _stencil_embed_face_value_x (point, u.x, 0);_stencil_val(u.x,0,0,0); _stencil_val(cs,0,0,0); _stencil_val(u.x,0 -1,0,0); _stencil_val(cs,0 -1,0,0);_stencil_val(cs,0,0,0); _stencil_val(cs,0 -1,0,0);                 }}end__stencil_is_face_x()
#line 191
_stencil_is_face_y(){
    {_stencil_val_a(uf.y,0,0,0); _stencil_val(fm.y,0,0,0); _stencil_val(fs.y,0,0,0); _stencil_val(fs.y,0,0,0); _stencil_embed_face_value_y (point, u.y, 0);_stencil_val(u.y,0,0,0); _stencil_val(cs,0,0,0); _stencil_val(u.y,0,0 -1,0); _stencil_val(cs,0,0 -1,0);_stencil_val(cs,0,0,0); _stencil_val(cs,0,0 -1,0);                 }}end__stencil_is_face_y()}end_foreach_face_stencil();
  
#line 191
if(!is_constant(fm.x)){{foreach_face_generic(){is_face_x(){
    val(uf.x,0,0,0) = val(fm.x,0,0,0)*(_attribute[u.x.i].third && val(fs.x,0,0,0) < 1. && val(fs.x,0,0,0) > 0. ? embed_face_value_x (point, u.x, 0) : ((val(u.x,0,0,0)*(1.5 + val(cs,0,0,0)) + val(u.x,0 -1,0,0)*(1.5 + val(cs,0 -1,0,0)))/ (val(cs,0,0,0) + val(cs,0 -1,0,0) + 3.)));}end_is_face_x()
#line 191
is_face_y(){
    val(uf.y,0,0,0) = val(fm.y,0,0,0)*(_attribute[u.y.i].third && val(fs.y,0,0,0) < 1. && val(fs.y,0,0,0) > 0. ? embed_face_value_y (point, u.y, 0) : ((val(u.y,0,0,0)*(1.5 + val(cs,0,0,0)) + val(u.y,0,0 -1,0)*(1.5 + val(cs,0,0 -1,0)))/ (val(cs,0,0,0) + val(cs,0,0 -1,0) + 3.)));}end_is_face_y()}end_foreach_face_generic();}}else {struct{double x,y;}_const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);
  {
#line 191
foreach_face_generic(){is_face_x(){
    val(uf.x,0,0,0) = _const_fm.x*(_attribute[u.x.i].third && val(fs.x,0,0,0) < 1. && val(fs.x,0,0,0) > 0. ? embed_face_value_x (point, u.x, 0) : ((val(u.x,0,0,0)*(1.5 + val(cs,0,0,0)) + val(u.x,0 -1,0,0)*(1.5 + val(cs,0 -1,0,0)))/ (val(cs,0,0,0) + val(cs,0 -1,0,0) + 3.)));}end_is_face_x()
#line 191
is_face_y(){
    val(uf.y,0,0,0) = _const_fm.y*(_attribute[u.y.i].third && val(fs.y,0,0,0) < 1. && val(fs.y,0,0,0) > 0. ? embed_face_value_y (point, u.y, 0) : ((val(u.y,0,0,0)*(1.5 + val(cs,0,0,0)) + val(u.y,0,0 -1,0)*(1.5 + val(cs,0,0 -1,0)))/ (val(cs,0,0,0) + val(cs,0,0 -1,0) + 3.)));}end_is_face_y()}end_foreach_face_generic();}}




  event ("properties");





  dtmax = DT;
  event ("stability");
}{end_tracing("init","/u/basilisk/src/navier-stokes/centered.h",0);return 0;}end_tracing("init","/u/basilisk/src/navier-stokes/centered.h",0);}
#line 214 "/u/basilisk/src/navier-stokes/centered.h"
static int set_dtmax_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i++);*ip=i;*tp=t;return ret;}      static int set_dtmax(const int i,const double t,Event *_ev){tracing("set_dtmax","/u/basilisk/src/navier-stokes/centered.h",0); dtmax = DT;{end_tracing("set_dtmax","/u/basilisk/src/navier-stokes/centered.h",0);return 0;}end_tracing("set_dtmax","/u/basilisk/src/navier-stokes/centered.h",0);}

static int stability_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i++);*ip=i;*tp=t;return ret;}      static int stability(const int i,const double t,Event *_ev){tracing("stability","/u/basilisk/src/navier-stokes/centered.h",0); {
  dt = dtnext (stokes ? dtmax : timestep (uf, dtmax));
}{end_tracing("stability","/u/basilisk/src/navier-stokes/centered.h",0);return 0;}end_tracing("stability","/u/basilisk/src/navier-stokes/centered.h",0);}







static int vof_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i++);*ip=i;*tp=t;return ret;}static int vof(const int i,const double t,Event *_ev){;return 0;}
static int tracer_advection_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i++);*ip=i;*tp=t;return ret;}static int tracer_advection(const int i,const double t,Event *_ev){;return 0;}
static int tracer_diffusion_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i++);*ip=i;*tp=t;return ret;}static int tracer_diffusion(const int i,const double t,Event *_ev){;return 0;}






static int properties_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i++);*ip=i;*tp=t;return ret;}static int properties(const int i,const double t,Event *_ev){;return 0;}
#line 247 "/u/basilisk/src/navier-stokes/centered.h"
void prediction()
{
  vector du;
   {
    scalar s = new_scalar("s");
    du.x = s;
  } 
#line 250
{
    scalar s = new_scalar("s");
    du.y = s;
  }

  if (_attribute[u.x.i].gradient)
    {
    
#line 256
foreach_stencil()
      { {

_stencil_val(fs.x,0,0,0);_stencil_val(fs.x,1,0,0);{
   {_stencil_val_a(du.x,0,0,0);  }

   
{_stencil_val_a(du.x,0,0,0);_stencil_val(u.x,-1,0,0); _stencil_val(u.x,0,0,0); _stencil_val(u.x,1,0,0);   }}

           
 
      
#line 264
} 
#line 257
{

_stencil_val(fs.y,0,0,0);_stencil_val(fs.y,0,1,0);{
   {_stencil_val_a(du.y,0,0,0);  }

   
{_stencil_val_a(du.y,0,0,0);_stencil_val(u.y,0,-1,0); _stencil_val(u.y,0,0,0); _stencil_val(u.y,0,1,0);   }}

           
 
      
#line 264
}}end_foreach_stencil();{
#line 256
foreach()
      { {

        if (!val(fs.x,0,0,0) || !val(fs.x,1,0,0))
   val(du.x,0,0,0) = 0.;
 else

   val(du.x,0,0,0) = _attribute[u.x.i].gradient (val(u.x,-1,0,0), val(u.x,0,0,0), val(u.x,1,0,0))/Delta;
      } 
#line 257
{

        if (!val(fs.y,0,0,0) || !val(fs.y,0,1,0))
   val(du.y,0,0,0) = 0.;
 else

   val(du.y,0,0,0) = _attribute[u.y.i].gradient (val(u.y,0,-1,0), val(u.y,0,0,0), val(u.y,0,1,0))/Delta;
      }}end_foreach();}}
  else
    {
    
#line 266
foreach_stencil()
      { {

_stencil_val(fs.x,0,0,0);_stencil_val(fs.x,1,0,0);{
   {_stencil_val_a(du.x,0,0,0);  }

   
{_stencil_val_a(du.x,0,0,0);_stencil_val(u.x,1,0,0); _stencil_val(u.x,-1,0,0);   }}

           
 
    
#line 274
} 
#line 267
{

_stencil_val(fs.y,0,0,0);_stencil_val(fs.y,0,1,0);{
   {_stencil_val_a(du.y,0,0,0);  }

   
{_stencil_val_a(du.y,0,0,0);_stencil_val(u.y,0,1,0); _stencil_val(u.y,0,-1,0);   }}

           
 
    
#line 274
}}end_foreach_stencil();{
#line 266
foreach()
      { {

        if (!val(fs.x,0,0,0) || !val(fs.x,1,0,0))
   val(du.x,0,0,0) = 0.;
 else

   val(du.x,0,0,0) = (val(u.x,1,0,0) - val(u.x,-1,0,0))/(2.*Delta);
    } 
#line 267
{

        if (!val(fs.y,0,0,0) || !val(fs.y,0,1,0))
   val(du.y,0,0,0) = 0.;
 else

   val(du.y,0,0,0) = (val(u.y,0,1,0) - val(u.y,0,-1,0))/(2.*Delta);
    }}end_foreach();}}

  trash (((vector[]){uf,{{-1},{-1}}}));
  foreach_face_stencil(){_stencil_is_face_x(){ {       
     _stencil_val(u.x,-1,0,0);_stencil_val(u.x,0,0,0);     
    
    _stencil_val_a(uf.x,0,0,0);_stencil_val(u.x, o_stencil,0,0);_stencil_val(g.x,0,0,0); _stencil_val(g.x,-1,0,0);_stencil_val(du.x,o_stencil,0,0);

_stencil_val(fm.y,o_stencil,0,0); _stencil_val(fm.y,o_stencil,1,0); {        
       _stencil_val(u.x,o_stencil,-1,0);_stencil_val(u.x, o_stencil,0,0);_stencil_val(u.x, o_stencil,0,0); _stencil_val(u.x,o_stencil,1,0);_stencil_val(u.y, o_stencil,0,0);
      _stencil_val_r(uf.x,0,0,0);_stencil_val(u.y,o_stencil,0,0);  
    }        

      







    
#line 293
_stencil_val_r(uf.x,0,0,0); _stencil_val(fm.x,0,0,0); 
  }}end__stencil_is_face_x()
#line 277
_stencil_is_face_y(){ {       
     _stencil_val(u.y,0,-1,0);_stencil_val(u.y,0,0,0);     
    
    _stencil_val_a(uf.y,0,0,0);_stencil_val(u.y,0, o_stencil,0);_stencil_val(g.y,0,0,0); _stencil_val(g.y,0,-1,0);_stencil_val(du.y,0,o_stencil,0);

_stencil_val(fm.x,0,o_stencil,0); _stencil_val(fm.x,1,o_stencil,0); {        
       _stencil_val(u.y,-1,o_stencil,0);_stencil_val(u.y,0, o_stencil,0);_stencil_val(u.y,0, o_stencil,0); _stencil_val(u.y,1,o_stencil,0);_stencil_val(u.x,0, o_stencil,0);
      _stencil_val_r(uf.y,0,0,0);_stencil_val(u.x,0,o_stencil,0);  
    }        

      







    
#line 293
_stencil_val_r(uf.y,0,0,0); _stencil_val(fm.y,0,0,0); 
  }}end__stencil_is_face_y()}end_foreach_face_stencil();
  
#line 277
if(!is_constant(fm.x)){{foreach_face_generic(){is_face_x(){ {
    double un = dt*(val(u.x,0,0,0) + val(u.x,-1,0,0))/(2.*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    val(uf.x,0,0,0) = val(u.x,i,0,0) + (val(g.x,0,0,0) + val(g.x,-1,0,0))*dt/4. + s*(1. - s*un)*val(du.x,i,0,0)*Delta/2.;

    if (val(fm.y,i,0,0) && val(fm.y,i,1,0)) {
      double fyy = val(u.y,i,0,0) < 0. ? val(u.x,i,1,0) - val(u.x,i,0,0) : val(u.x,i,0,0) - val(u.x,i,-1,0);
      val(uf.x,0,0,0) -= dt*val(u.y,i,0,0)*fyy/(2.*Delta);
    }







    val(uf.x,0,0,0) *= val(fm.x,0,0,0);
  }}end_is_face_x()
#line 277
is_face_y(){ {
    double un = dt*(val(u.y,0,0,0) + val(u.y,0,-1,0))/(2.*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    val(uf.y,0,0,0) = val(u.y,0,i,0) + (val(g.y,0,0,0) + val(g.y,0,-1,0))*dt/4. + s*(1. - s*un)*val(du.y,0,i,0)*Delta/2.;

    if (val(fm.x,0,i,0) && val(fm.x,1,i,0)) {
      double fyy = val(u.x,0,i,0) < 0. ? val(u.y,1,i,0) - val(u.y,0,i,0) : val(u.y,0,i,0) - val(u.y,-1,i,0);
      val(uf.y,0,0,0) -= dt*val(u.x,0,i,0)*fyy/(2.*Delta);
    }







    val(uf.y,0,0,0) *= val(fm.y,0,0,0);
  }}end_is_face_y()}end_foreach_face_generic();}}else {struct{double x,y;}_const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);
  {
#line 277
foreach_face_generic(){is_face_x(){ {
    double un = dt*(val(u.x,0,0,0) + val(u.x,-1,0,0))/(2.*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    val(uf.x,0,0,0) = val(u.x,i,0,0) + (val(g.x,0,0,0) + val(g.x,-1,0,0))*dt/4. + s*(1. - s*un)*val(du.x,i,0,0)*Delta/2.;

    if (_const_fm.y && _const_fm.y) {
      double fyy = val(u.y,i,0,0) < 0. ? val(u.x,i,1,0) - val(u.x,i,0,0) : val(u.x,i,0,0) - val(u.x,i,-1,0);
      val(uf.x,0,0,0) -= dt*val(u.y,i,0,0)*fyy/(2.*Delta);
    }







    val(uf.x,0,0,0) *= _const_fm.x;
  }}end_is_face_x()
#line 277
is_face_y(){ {
    double un = dt*(val(u.y,0,0,0) + val(u.y,0,-1,0))/(2.*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    val(uf.y,0,0,0) = val(u.y,0,i,0) + (val(g.y,0,0,0) + val(g.y,0,-1,0))*dt/4. + s*(1. - s*un)*val(du.y,0,i,0)*Delta/2.;

    if (_const_fm.x && _const_fm.x) {
      double fyy = val(u.x,0,i,0) < 0. ? val(u.y,1,i,0) - val(u.y,0,i,0) : val(u.y,0,i,0) - val(u.y,-1,i,0);
      val(uf.y,0,0,0) -= dt*val(u.x,0,i,0)*fyy/(2.*Delta);
    }







    val(uf.y,0,0,0) *= _const_fm.y;
  }}end_is_face_y()}end_foreach_face_generic();}}

  delete ((scalar *)((vector[]){du,{{-1},{-1}}}));
}
#line 308 "/u/basilisk/src/navier-stokes/centered.h"
static int advection_term_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i++);*ip=i;*tp=t;return ret;}      static int advection_term(const int i,const double t,Event *_ev){tracing("advection_term","/u/basilisk/src/navier-stokes/centered.h",0);
{
  if (!stokes) {
    prediction();
    mgpf = project ((struct Project){uf, pf, alpha, dt/2., mgpf.nrelax});
    advection ((struct Advection){(scalar *)((vector[]){u,{{-1},{-1}}}), uf, dt, (scalar *)((vector[]){g,{{-1},{-1}}})});
  }
}{end_tracing("advection_term","/u/basilisk/src/navier-stokes/centered.h",0);return 0;}end_tracing("advection_term","/u/basilisk/src/navier-stokes/centered.h",0);}







static void correction (double dt)
{
  foreach_stencil()
    {
      {_stencil_val_r(u.x,0,0,0);_stencil_val(g.x,0,0,0);  }
      
#line 327
{_stencil_val_r(u.y,0,0,0);_stencil_val(g.y,0,0,0);  }}end_foreach_stencil();
  {
#line 325
foreach()
    {
      val(u.x,0,0,0) += dt*val(g.x,0,0,0);
      
#line 327
val(u.y,0,0,0) += dt*val(g.y,0,0,0);}end_foreach();}
}
#line 337 "/u/basilisk/src/navier-stokes/centered.h"
static int viscous_term_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i++);*ip=i;*tp=t;return ret;}      static int viscous_term(const int i,const double t,Event *_ev){tracing("viscous_term","/u/basilisk/src/navier-stokes/centered.h",0);
{
  if (constant(mu.x) != 0.) {
    correction (dt);
    mgu = viscosity ((struct Viscosity){u, mu, rho, dt, mgu.nrelax});
    correction (-dt);
  }




  if (!is_constant(a.x)) {
    vector af = a;
    trash (((vector[]){af,{{-1},{-1}}}));
    foreach_face_stencil(){_stencil_is_face_x(){
      {_stencil_val_a(af.x,0,0,0);  }}end__stencil_is_face_x()
#line 351
_stencil_is_face_y(){
      {_stencil_val_a(af.y,0,0,0);  }}end__stencil_is_face_y()}end_foreach_face_stencil();
    {
#line 351
foreach_face_generic(){is_face_x(){
      val(af.x,0,0,0) = 0.;}end_is_face_x()
#line 351
is_face_y(){
      val(af.y,0,0,0) = 0.;}end_is_face_y()}end_foreach_face_generic();}
  }
}{end_tracing("viscous_term","/u/basilisk/src/navier-stokes/centered.h",0);return 0;}end_tracing("viscous_term","/u/basilisk/src/navier-stokes/centered.h",0);}
#line 373 "/u/basilisk/src/navier-stokes/centered.h"
static int acceleration_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i++);*ip=i;*tp=t;return ret;}      static int acceleration(const int i,const double t,Event *_ev){tracing("acceleration","/u/basilisk/src/navier-stokes/centered.h",0);
{
  trash (((vector[]){uf,{{-1},{-1}}}));
  foreach_face_stencil(){_stencil_is_face_x(){
    {_stencil_val_a(uf.x,0,0,0); _stencil_val(fm.x,0,0,0); _stencil_val(fs.x,0,0,0); _stencil_val(fs.x,0,0,0); _stencil_embed_face_value_x (point, u.x, 0);_stencil_val(u.x,0,0,0); _stencil_val(cs,0,0,0); _stencil_val(u.x,0 -1,0,0); _stencil_val(cs,0 -1,0,0);_stencil_val(cs,0,0,0); _stencil_val(cs,0 -1,0,0);_stencil_val(a.x,0,0,0);                   }}end__stencil_is_face_x()
#line 376
_stencil_is_face_y(){
    {_stencil_val_a(uf.y,0,0,0); _stencil_val(fm.y,0,0,0); _stencil_val(fs.y,0,0,0); _stencil_val(fs.y,0,0,0); _stencil_embed_face_value_y (point, u.y, 0);_stencil_val(u.y,0,0,0); _stencil_val(cs,0,0,0); _stencil_val(u.y,0,0 -1,0); _stencil_val(cs,0,0 -1,0);_stencil_val(cs,0,0,0); _stencil_val(cs,0,0 -1,0);_stencil_val(a.y,0,0,0);                   }}end__stencil_is_face_y()}end_foreach_face_stencil();
  
#line 376
if(!is_constant(fm.x) && !is_constant(a.x)){{foreach_face_generic(){is_face_x(){
    val(uf.x,0,0,0) = val(fm.x,0,0,0)*((_attribute[u.x.i].third && val(fs.x,0,0,0) < 1. && val(fs.x,0,0,0) > 0. ? embed_face_value_x (point, u.x, 0) : ((val(u.x,0,0,0)*(1.5 + val(cs,0,0,0)) + val(u.x,0 -1,0,0)*(1.5 + val(cs,0 -1,0,0)))/ (val(cs,0,0,0) + val(cs,0 -1,0,0) + 3.))) + dt*val(a.x,0,0,0));}end_is_face_x()
#line 376
is_face_y(){
    val(uf.y,0,0,0) = val(fm.y,0,0,0)*((_attribute[u.y.i].third && val(fs.y,0,0,0) < 1. && val(fs.y,0,0,0) > 0. ? embed_face_value_y (point, u.y, 0) : ((val(u.y,0,0,0)*(1.5 + val(cs,0,0,0)) + val(u.y,0,0 -1,0)*(1.5 + val(cs,0,0 -1,0)))/ (val(cs,0,0,0) + val(cs,0,0 -1,0) + 3.))) + dt*val(a.y,0,0,0));}end_is_face_y()}end_foreach_face_generic();}}else if(is_constant(fm.x) && !is_constant(a.x)){struct{double x,y;}_const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);
  {
#line 376
foreach_face_generic(){is_face_x(){
    val(uf.x,0,0,0) = _const_fm.x*((_attribute[u.x.i].third && val(fs.x,0,0,0) < 1. && val(fs.x,0,0,0) > 0. ? embed_face_value_x (point, u.x, 0) : ((val(u.x,0,0,0)*(1.5 + val(cs,0,0,0)) + val(u.x,0 -1,0,0)*(1.5 + val(cs,0 -1,0,0)))/ (val(cs,0,0,0) + val(cs,0 -1,0,0) + 3.))) + dt*val(a.x,0,0,0));}end_is_face_x()
#line 376
is_face_y(){
    val(uf.y,0,0,0) = _const_fm.y*((_attribute[u.y.i].third && val(fs.y,0,0,0) < 1. && val(fs.y,0,0,0) > 0. ? embed_face_value_y (point, u.y, 0) : ((val(u.y,0,0,0)*(1.5 + val(cs,0,0,0)) + val(u.y,0,0 -1,0)*(1.5 + val(cs,0,0 -1,0)))/ (val(cs,0,0,0) + val(cs,0,0 -1,0) + 3.))) + dt*val(a.y,0,0,0));}end_is_face_y()}end_foreach_face_generic();}}else if(!is_constant(fm.x) && is_constant(a.x)){struct{double x,y;}_const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);
  {
#line 376
foreach_face_generic(){is_face_x(){
    val(uf.x,0,0,0) = val(fm.x,0,0,0)*((_attribute[u.x.i].third && val(fs.x,0,0,0) < 1. && val(fs.x,0,0,0) > 0. ? embed_face_value_x (point, u.x, 0) : ((val(u.x,0,0,0)*(1.5 + val(cs,0,0,0)) + val(u.x,0 -1,0,0)*(1.5 + val(cs,0 -1,0,0)))/ (val(cs,0,0,0) + val(cs,0 -1,0,0) + 3.))) + dt*_const_a.x);}end_is_face_x()
#line 376
is_face_y(){
    val(uf.y,0,0,0) = val(fm.y,0,0,0)*((_attribute[u.y.i].third && val(fs.y,0,0,0) < 1. && val(fs.y,0,0,0) > 0. ? embed_face_value_y (point, u.y, 0) : ((val(u.y,0,0,0)*(1.5 + val(cs,0,0,0)) + val(u.y,0,0 -1,0)*(1.5 + val(cs,0,0 -1,0)))/ (val(cs,0,0,0) + val(cs,0,0 -1,0) + 3.))) + dt*_const_a.y);}end_is_face_y()}end_foreach_face_generic();}}else {struct{double x,y;}_const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);struct{double x,y;}_const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);
  {
#line 376
foreach_face_generic(){is_face_x(){
    val(uf.x,0,0,0) = _const_fm.x*((_attribute[u.x.i].third && val(fs.x,0,0,0) < 1. && val(fs.x,0,0,0) > 0. ? embed_face_value_x (point, u.x, 0) : ((val(u.x,0,0,0)*(1.5 + val(cs,0,0,0)) + val(u.x,0 -1,0,0)*(1.5 + val(cs,0 -1,0,0)))/ (val(cs,0,0,0) + val(cs,0 -1,0,0) + 3.))) + dt*_const_a.x);}end_is_face_x()
#line 376
is_face_y(){
    val(uf.y,0,0,0) = _const_fm.y*((_attribute[u.y.i].third && val(fs.y,0,0,0) < 1. && val(fs.y,0,0,0) > 0. ? embed_face_value_y (point, u.y, 0) : ((val(u.y,0,0,0)*(1.5 + val(cs,0,0,0)) + val(u.y,0,0 -1,0)*(1.5 + val(cs,0,0 -1,0)))/ (val(cs,0,0,0) + val(cs,0,0 -1,0) + 3.))) + dt*_const_a.y);}end_is_face_y()}end_foreach_face_generic();}}
}{end_tracing("acceleration","/u/basilisk/src/navier-stokes/centered.h",0);return 0;}end_tracing("acceleration","/u/basilisk/src/navier-stokes/centered.h",0);}
#line 387 "/u/basilisk/src/navier-stokes/centered.h"
void centered_gradient (scalar p, vector g)
{





  vector  gf=new_face_vector("gf");
  foreach_face_stencil(){_stencil_is_face_x(){
    {_stencil_val_a(gf.x,0,0,0); _stencil_val(fm.x,0,0,0);_stencil_val(a.x,0,0,0); _stencil_val(alpha.x,0,0,0);_stencil_val(p,0,0,0); _stencil_val(p,-1,0,0);   }}end__stencil_is_face_x()
#line 395
_stencil_is_face_y(){
    {_stencil_val_a(gf.y,0,0,0); _stencil_val(fm.y,0,0,0);_stencil_val(a.y,0,0,0); _stencil_val(alpha.y,0,0,0);_stencil_val(p,0,0,0); _stencil_val(p,0,-1,0);   }}end__stencil_is_face_y()}end_foreach_face_stencil();
  
#line 395
if(!is_constant(fm.x) && !is_constant(a.x) && !is_constant(alpha.x)){{foreach_face_generic(){is_face_x(){
    val(gf.x,0,0,0) = val(fm.x,0,0,0)*val(a.x,0,0,0) - val(alpha.x,0,0,0)*(val(p,0,0,0) - val(p,-1,0,0))/Delta;}end_is_face_x()
#line 395
is_face_y(){
    val(gf.y,0,0,0) = val(fm.y,0,0,0)*val(a.y,0,0,0) - val(alpha.y,0,0,0)*(val(p,0,0,0) - val(p,0,-1,0))/Delta;}end_is_face_y()}end_foreach_face_generic();}}else if(is_constant(fm.x) && !is_constant(a.x) && !is_constant(alpha.x)){struct{double x,y;}_const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);
  {
#line 395
foreach_face_generic(){is_face_x(){
    val(gf.x,0,0,0) = _const_fm.x*val(a.x,0,0,0) - val(alpha.x,0,0,0)*(val(p,0,0,0) - val(p,-1,0,0))/Delta;}end_is_face_x()
#line 395
is_face_y(){
    val(gf.y,0,0,0) = _const_fm.y*val(a.y,0,0,0) - val(alpha.y,0,0,0)*(val(p,0,0,0) - val(p,0,-1,0))/Delta;}end_is_face_y()}end_foreach_face_generic();}}else if(!is_constant(fm.x) && is_constant(a.x) && !is_constant(alpha.x)){struct{double x,y;}_const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);
  {
#line 395
foreach_face_generic(){is_face_x(){
    val(gf.x,0,0,0) = val(fm.x,0,0,0)*_const_a.x - val(alpha.x,0,0,0)*(val(p,0,0,0) - val(p,-1,0,0))/Delta;}end_is_face_x()
#line 395
is_face_y(){
    val(gf.y,0,0,0) = val(fm.y,0,0,0)*_const_a.y - val(alpha.y,0,0,0)*(val(p,0,0,0) - val(p,0,-1,0))/Delta;}end_is_face_y()}end_foreach_face_generic();}}else if(is_constant(fm.x) && is_constant(a.x) && !is_constant(alpha.x)){struct{double x,y;}_const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);struct{double x,y;}_const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);
  {
#line 395
foreach_face_generic(){is_face_x(){
    val(gf.x,0,0,0) = _const_fm.x*_const_a.x - val(alpha.x,0,0,0)*(val(p,0,0,0) - val(p,-1,0,0))/Delta;}end_is_face_x()
#line 395
is_face_y(){
    val(gf.y,0,0,0) = _const_fm.y*_const_a.y - val(alpha.y,0,0,0)*(val(p,0,0,0) - val(p,0,-1,0))/Delta;}end_is_face_y()}end_foreach_face_generic();}}else if(!is_constant(fm.x) && !is_constant(a.x) && is_constant(alpha.x)){struct{double x,y;}_const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);
  {
#line 395
foreach_face_generic(){is_face_x(){
    val(gf.x,0,0,0) = val(fm.x,0,0,0)*val(a.x,0,0,0) - _const_alpha.x*(val(p,0,0,0) - val(p,-1,0,0))/Delta;}end_is_face_x()
#line 395
is_face_y(){
    val(gf.y,0,0,0) = val(fm.y,0,0,0)*val(a.y,0,0,0) - _const_alpha.y*(val(p,0,0,0) - val(p,0,-1,0))/Delta;}end_is_face_y()}end_foreach_face_generic();}}else if(is_constant(fm.x) && !is_constant(a.x) && is_constant(alpha.x)){struct{double x,y;}_const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);struct{double x,y;}_const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);
  {
#line 395
foreach_face_generic(){is_face_x(){
    val(gf.x,0,0,0) = _const_fm.x*val(a.x,0,0,0) - _const_alpha.x*(val(p,0,0,0) - val(p,-1,0,0))/Delta;}end_is_face_x()
#line 395
is_face_y(){
    val(gf.y,0,0,0) = _const_fm.y*val(a.y,0,0,0) - _const_alpha.y*(val(p,0,0,0) - val(p,0,-1,0))/Delta;}end_is_face_y()}end_foreach_face_generic();}}else if(!is_constant(fm.x) && is_constant(a.x) && is_constant(alpha.x)){struct{double x,y;}_const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);struct{double x,y;}_const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);
  {
#line 395
foreach_face_generic(){is_face_x(){
    val(gf.x,0,0,0) = val(fm.x,0,0,0)*_const_a.x - _const_alpha.x*(val(p,0,0,0) - val(p,-1,0,0))/Delta;}end_is_face_x()
#line 395
is_face_y(){
    val(gf.y,0,0,0) = val(fm.y,0,0,0)*_const_a.y - _const_alpha.y*(val(p,0,0,0) - val(p,0,-1,0))/Delta;}end_is_face_y()}end_foreach_face_generic();}}else {struct{double x,y;}_const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);struct{double x,y;}_const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);struct{double x,y;}_const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);
  {
#line 395
foreach_face_generic(){is_face_x(){
    val(gf.x,0,0,0) = _const_fm.x*_const_a.x - _const_alpha.x*(val(p,0,0,0) - val(p,-1,0,0))/Delta;}end_is_face_x()
#line 395
is_face_y(){
    val(gf.y,0,0,0) = _const_fm.y*_const_a.y - _const_alpha.y*(val(p,0,0,0) - val(p,0,-1,0))/Delta;}end_is_face_y()}end_foreach_face_generic();}}





  trash (((vector[]){g,{{-1},{-1}}}));
  foreach_stencil()
    {
      {_stencil_val_a(g.x,0,0,0);_stencil_val(gf.x,0,0,0); _stencil_val(gf.x,1,0,0);_stencil_val(fm.x,0,0,0); _stencil_val(fm.x,1,0,0);      }
      
#line 405
{_stencil_val_a(g.y,0,0,0);_stencil_val(gf.y,0,0,0); _stencil_val(gf.y,0,1,0);_stencil_val(fm.y,0,0,0); _stencil_val(fm.y,0,1,0);      }}end_foreach_stencil();
  
#line 403
if(!is_constant(fm.x)){{foreach()
    {
      val(g.x,0,0,0) = (val(gf.x,0,0,0) + val(gf.x,1,0,0))/(val(fm.x,0,0,0) + val(fm.x,1,0,0) + 1e-30);
      
#line 405
val(g.y,0,0,0) = (val(gf.y,0,0,0) + val(gf.y,0,1,0))/(val(fm.y,0,0,0) + val(fm.y,0,1,0) + 1e-30);}end_foreach();}}else {struct{double x,y;}_const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);
  {
#line 403
foreach()
    {
      val(g.x,0,0,0) = (val(gf.x,0,0,0) + val(gf.x,1,0,0))/(_const_fm.x + _const_fm.x + 1e-30);
      
#line 405
val(g.y,0,0,0) = (val(gf.y,0,0,0) + val(gf.y,0,1,0))/(_const_fm.y + _const_fm.y + 1e-30);}end_foreach();}}delete((scalar*)((vector[]){gf,{{-1},{-1}}}));
}






static int projection_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i++);*ip=i;*tp=t;return ret;}      static int projection(const int i,const double t,Event *_ev){tracing("projection","/u/basilisk/src/navier-stokes/centered.h",0);
{
  mgp = project ((struct Project){uf, p, alpha, dt, mgp.nrelax});
  centered_gradient (p, g);




  correction (dt);
}{end_tracing("projection","/u/basilisk/src/navier-stokes/centered.h",0);return 0;}end_tracing("projection","/u/basilisk/src/navier-stokes/centered.h",0);}





static int end_timestep_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i++);*ip=i;*tp=t;return ret;}static int end_timestep(const int i,const double t,Event *_ev){;return 0;}
#line 438 "/u/basilisk/src/navier-stokes/centered.h"
static int adapt_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i++);*ip=i;*tp=t;return ret;}      static int adapt(const int i,const double t,Event *_ev){tracing("adapt","/u/basilisk/src/navier-stokes/centered.h",0); {

  fractions_cleanup ((struct Cleanup){cs, fs});
  foreach_face_stencil(){_stencil_is_face_x(){
    {_stencil_val(uf.x,0,0,0);_stencil_val(fs.x,0,0,0);
      {_stencil_val_a(uf.x,0,0,0);  }   }}end__stencil_is_face_x()
#line 441
_stencil_is_face_y(){
    {_stencil_val(uf.y,0,0,0);_stencil_val(fs.y,0,0,0);
      {_stencil_val_a(uf.y,0,0,0);  }   }}end__stencil_is_face_y()}end_foreach_face_stencil();
  {
#line 441
foreach_face_generic(){is_face_x(){
    if (val(uf.x,0,0,0) && !val(fs.x,0,0,0))
      val(uf.x,0,0,0) = 0.;}end_is_face_x()
#line 441
is_face_y(){
    if (val(uf.y,0,0,0) && !val(fs.y,0,0,0))
      val(uf.y,0,0,0) = 0.;}end_is_face_y()}end_foreach_face_generic();}

  event ("properties");
}{end_tracing("adapt","/u/basilisk/src/navier-stokes/centered.h",0);return 0;}end_tracing("adapt","/u/basilisk/src/navier-stokes/centered.h",0);}
#line 5 "cylinder.c"
#line 1 "output_htg.h"
#line 1 "./output_htg.h"
void output_htg(scalar *list, vector *vlist, const char *path);




void output_htg_data(scalar *list, vector *vlist, FILE *fp);
#line 38 "./output_htg.h"
void output_htg(scalar *list, vector *vlist, const char *path) {
  FILE *fp;
  fp = fopen(path, "w");
  if (!fp) {
    printf("output_htg.h : %s could not be opened\n Does the Folder exist?\n",
           path);
    exit(1);
  }

  output_htg_data((scalar *)list, (vector *)vlist, fp);
  fclose(fp);
}
#line 592 "./output_htg.h"
void output_htg_data(scalar *list, vector *vlist, FILE *fp) {
  unsigned int vertices_local = 0;
  unsigned int descBits_local = 0;
  unsigned int vertices_local_pL[grid->maxdepth + 1];

  for (int lvl = 0; lvl < grid->maxdepth + 1; ++lvl) {
    vertices_local_pL[lvl] = 0;
    
#if _OPENMP
  #undef OMP
  #define OMP(x)
#endif
{
#line 599
foreach_level(lvl) if (is_local(cell)) vertices_local_pL[lvl]++;end_foreach_level();}
#if _OPENMP
  #undef OMP
  #define OMP(x) _Pragma(#x)
#endif


    
#line 601
vertices_local += vertices_local_pL[lvl];
  }

  descBits_local = vertices_local - vertices_local_pL[grid->maxdepth];

  double min_val[list_len(list)];
  double max_val[list_len(list)];
  {
    int i = 0;
    {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {
      stats stat = statsf(s);
      min_val[i] = stat.min;
      max_val[i] = stat.max;
      i++;
    }}}
  }
  double min_val_v[vectors_len(vlist)];
  double max_val_v[vectors_len(vlist)];
  {
    int i = 0;
    {vector*_i=(vector*)( vlist);if(_i)for(vector v=*_i;(&v)->x.i>=0;v=*++_i){ {
      min_val_v[i] = HUGE;
      max_val_v[i] = -HUGE;
       {
        stats stat = statsf(v.x);
        min_val_v[i] = min(stat.min, min_val_v[i]);
        max_val_v[i] = max(stat.max, max_val_v[i]);
      } 
#line 624
{
        stats stat = statsf(v.y);
        min_val_v[i] = min(stat.min, min_val_v[i]);
        max_val_v[i] = max(stat.max, max_val_v[i]);
      }
      i++;
    }}}
  }
  int maj_v = 1, min_v = 0;

  fprintf(fp, "<VTKFile %s version=\"%i.%i\" %s %s>\n",
          "type=\"HyperTreeGrid\"", maj_v, min_v,
          "byte_order=\"LittleEndian\" ", "header_type=\"UInt32\"");


  fprintf(fp,
          "\t<HyperTreeGrid BranchFactor=\"2\" TransposedRootIndexing=\"0\" "
          "Dimensions=\"%d %d %d\">\n",
          2, 2, 1);






  fprintf(fp, "\t\t<Grid>\n");

  fprintf(fp,
          "\t\t\t<DataArray type=\"Float64\" Name=\"XCoordinates\" "
          "NumberOfTuples=\"2\" format=\"ascii\" RangeMin=\"%g\" "
          "RangeMax=\"%g\">\n",
          Y0, Y0 + L0);
  fprintf(fp, "\t\t\t\t%g %g", Y0, Y0 + L0);
  fprintf(fp, "\n\t\t\t</DataArray>\n");
  fprintf(fp,
          "\t\t\t<DataArray type=\"Float64\" Name=\"YCoordinates\" "
          "NumberOfTuples=\"2\" format=\"ascii\" RangeMin=\"%g\" "
          "RangeMax=\"%g\">\n",
          X0, X0 + L0);
  fprintf(fp, "\t\t\t\t%g %g", X0, X0 + L0);
  fprintf(fp, "\n\t\t\t</DataArray>\n");
  fprintf(fp,
          "\t\t\t<DataArray type=\"Float64\" Name=\"ZCoordinates\" "
          "NumberOfTuples=\"2\" format=\"ascii\" RangeMin=\"%g\" "
          "RangeMax=\"%g\">\n",
          Z0, Z0 + L0);
  fprintf(fp, "\t\t\t\t%g %g", 0., 0.);
#line 693 "./output_htg.h"
  fprintf(fp, "\n\t\t\t</DataArray>\n");
  fprintf(fp, "\t\t</Grid>\n");
  fprintf(fp, "\t\t<Trees>\n");

  unsigned int byte_offset = 0;
  fprintf(fp,
          "\t\t\t<Tree Index=\"0\" NumberOfLevels=\"%d\" "
          "NumberOfVertices=\"%u\">\n",
          grid->maxdepth + 1, vertices_local);

  fprintf(fp,
          "\t\t\t\t<DataArray type=\"Bit\" Name=\"Descriptor\" "
          "NumberOfTuples=\"%u\" format=\"appended\" RangeMin=\"0\" "
          "RangeMax=\"1\" offset=\"%u\"/>\n",
          descBits_local, byte_offset);
  byte_offset += (descBits_local / 8 + 1) * sizeof(uint8_t) + sizeof(uint32_t);

  fprintf(fp,
          "\t\t\t\t<DataArray type=\"Int64\" Name=\"NbVerticesByLevel\" "
          "NumberOfTuples=\"%d\" format=\"ascii\" RangeMin=\"1\" "
          "RangeMax=\"%u\" >\n\t\t\t\t\t",
          grid->maxdepth + 1, vertices_local_pL[grid->maxdepth]);

  for (int lvl = 0; lvl <= grid->maxdepth; lvl++) {
    fprintf(fp, "%u ", vertices_local_pL[lvl]);
  }
  fprintf(fp, "\n\t\t\t\t</DataArray>\n");
  fprintf(fp, "\t\t\t\t<CellData>\n");
  {
    int i = 0;
    {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {
      fprintf(fp,
              "\t\t\t\t\t<DataArray type=\"Float32\" Name=\"%s\" "
              "NumberOfTuples=\"%u\" format=\"appended\" RangeMin=\"%g\" "
              "RangeMax=\"%g\" offset=\"%u\"/>\n",
              _attribute[s.i].name, vertices_local, min_val[i], max_val[i], byte_offset);
      byte_offset += vertices_local * sizeof(float) + sizeof(uint32_t);
      i++;
    }}}
  }
  {
    int i = 0;
    {vector*_i=(vector*)( vlist);if(_i)for(vector v=*_i;(&v)->x.i>=0;v=*++_i){ {
      char *vname = strtok(_attribute[v.x.i].name, ".");
      fprintf(fp,
              "\t\t\t\t\t<DataArray type=\"Float32\" NumberOfComponents=\"%i\" "
              "Name=\"%s\" NumberOfTuples=\"%u\" format=\"appended\"  "
              "RangeMin=\"%g\" RangeMax=\"%g\" offset=\"%u\"/>\n",
              3, vname, vertices_local, min_val_v[i], max_val_v[i],
              byte_offset);
      byte_offset += vertices_local * 3 * sizeof(float) + sizeof(uint32_t);
      i++;
    }}}
  }
  fprintf(fp, "\t\t\t\t</CellData>\n");
  fprintf(fp, "\t\t\t</Tree>\n\t\t</Trees>\n");
  fprintf(fp, "\t</HyperTreeGrid>\n\t<AppendedData encoding=\"raw\">\n_");

  int cell_size;

  cell_size = sizeof(uint8_t);

  int vertices_local_corr = ((descBits_local / 8) + 1) * 8;

  uint32_t prepend_size = vertices_local_corr;
  fwrite(&prepend_size, sizeof(uint32_t), 1, fp);
  uint8_t *write_cache = (uint8_t *)calloc(vertices_local_corr, cell_size);
  long index = 1;
  for (int lvl = 0; lvl < grid->maxdepth; ++lvl) {
    
#if _OPENMP
  #undef OMP
  #define OMP(x)
#endif
{
#line 762
foreach_level(lvl) if (is_local(cell)) {
      if (is_leaf(cell)) {
        write_cache[index++] = 0;
      } else {
        write_cache[index++] = 1;
      }
    }end_foreach_level();}
#if _OPENMP
  #undef OMP
  #define OMP(x) _Pragma(#x)
#endif

  
#line 769
}

  for (int i = 0; i < vertices_local_corr / 8; ++i) {
    for (int j = 0; j < 8; ++j) {
      write_cache[i] |= write_cache[(8 * i + j) + 1] << (7 - j);
      if ((j + 1) == 8) {
        write_cache[i + 1] = 0;
      }
    }
  }

  fwrite(&write_cache[0], cell_size, vertices_local_corr / 8, fp);
  pfree(write_cache,__func__,__FILE__,0);
  write_cache = NULL;
  {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {
    cell_size = sizeof(float);

    uint32_t prepend_size = vertices_local * cell_size;
    fwrite(&prepend_size, sizeof(uint32_t), 1, fp);

    for (int lvl = 0; lvl < grid->maxdepth + 1; ++lvl) {

      float *write_cache = pmalloc(vertices_local_pL[lvl] * cell_size,__func__,__FILE__,0);
      long index = 0;

      
#if _OPENMP
  #undef OMP
  #define OMP(x)
#endif
{
#line 794
foreach_level(lvl) if (is_local(cell)) write_cache[index++] =
          val(s,0,0,0);end_foreach_level();}
#if _OPENMP
  #undef OMP
  #define OMP(x) _Pragma(#x)
#endif


      
#line 797
fwrite(&write_cache[0], cell_size, vertices_local_pL[lvl], fp);
      pfree(write_cache,__func__,__FILE__,0);
      write_cache = NULL;
    }
  }}}

  {vector*_i=(vector*)( vlist);if(_i)for(vector v=*_i;(&v)->x.i>=0;v=*++_i){ {
    cell_size = 3 * sizeof(float);

    uint32_t prepend_size = vertices_local * cell_size;
    fwrite(&prepend_size, sizeof(uint32_t), 1, fp);

    for (int lvl = 0; lvl <= grid->maxdepth; ++lvl) {

      float *write_cache = pmalloc(vertices_local_pL[lvl] * cell_size,__func__,__FILE__,0);
      long index = 0;
      
#if _OPENMP
  #undef OMP
  #define OMP(x)
#endif
{
#line 813
foreach_level(lvl) if (is_local(cell)) {


        write_cache[index] = val(v.y,0,0,0);
        write_cache[index + 1] = val(v.x,0,0,0);
        write_cache[index + 2] = 0.;
        index += 3;







      }end_foreach_level();}
#if _OPENMP
  #undef OMP
  #define OMP(x) _Pragma(#x)
#endif

      
#line 828
fwrite(&write_cache[0], cell_size, vertices_local_pL[lvl], fp);
      pfree(write_cache,__func__,__FILE__,0);
      write_cache = NULL;
    }
  }}}

  fprintf(fp, "ENDBINARY\n\t</AppendedData>\n</VTKFile>\n");
  fflush(fp);
}
#line 6 "cylinder.c"
static const double diameter = 0.125;
static const int minlevel = 5;
static double reynolds, tend;
static int level, period, Image, Surface;
static double _boundary4(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return dirichlet(1.);}}static double _boundary4_homogeneous(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return dirichlet_homogeneous();}}
static double _boundary5(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann(0.);}}static double _boundary5_homogeneous(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann_homogeneous();}}
static double _boundary6(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann(0.);}}static double _boundary6_homogeneous(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann_homogeneous();}}
static double _boundary7(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann(0.);}}static double _boundary7_homogeneous(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return neumann_homogeneous();}}
static double _boundary8(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return dirichlet(0.);}}static double _boundary8_homogeneous(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return dirichlet_homogeneous();}}
static double _boundary9(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return dirichlet(0.);}}static double _boundary9_homogeneous(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return dirichlet_homogeneous();}}
static double _boundary10(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return fabs(y) > 0.45 ? neumann(0.) : dirichlet(0.);}}static double _boundary10_homogeneous(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return fabs(y) > 0.45 ? neumann_homogeneous() : dirichlet_homogeneous();}}
static double _boundary11(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return fabs(y) > 0.45 ? neumann(0.) : dirichlet(0.);}}static double _boundary11_homogeneous(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return fabs(y) > 0.45 ? neumann_homogeneous() : dirichlet_homogeneous();}}
vector  muv={{11},{12}};
int main(int argc, char **argv) {_init_solver();
  char *end;
  int ReynoldsFlag;
  int LevelFlag;
  int PeriodFlag;
  int TendFlag;
  ReynoldsFlag = 0;
  LevelFlag = 0;
  PeriodFlag = 0;
  Image = 0;
  TendFlag = 0;
  while (*++argv != NULL && argv[0][0] == '-')
    switch (argv[0][1]) {
    case 'h':
      fprintf(
          ferr,
          "Usage: cylinder [-h] [-i] [-s] [-z <number of cells>] -r <Reynolds "
          "number> -l <resolution level> -p <dump period>\n"
          "Options:\n"
          "  -h     Display this help message\n"
          "  -i     Enable PPM image dumping\n"
          "  -r <Reynolds number>     the Reynolds number (a decimal number)\n"
          "  -l <resolution level>    the resolution level (positive integer)\n"
          "  -p <dump period>         the dump period (positive integer)\n"
          "  -e <end time>            end time of the simulation (decimal "
          "number)\n"
          "\n"
          "Example usage:\n"
          "  ./cylinder -i -r 100 -l 10 -p 100\n"
          "  ./cylinder -r 100 -l 10 -p 100\n");
      exit(1);
    case 'r':
      argv++;
      if (*argv == NULL) {
        fprintf(ferr, "cylinder: error:  -r needs an argument\n");
        exit(1);
      }
      reynolds = strtod(*argv, &end);
      if (*end != '\0') {
        fprintf(ferr, "cylinder: error: '%s' is not a number\n", *argv);
        exit(1);
      }
      ReynoldsFlag = 1;
      break;
    case 'l':
      argv++;
      if (*argv == NULL) {
        fprintf(ferr, "cylinder: error: -l needs an argument\n");
        exit(1);
      }
      level = strtol(*argv, &end, 10);
      if (*end != '\0' || level <= 0) {
        fprintf(ferr, "cylinder: error: '%s' is not a positive integer\n",
                *argv);
        exit(1);
      }
      LevelFlag = 1;
      break;
    case 'p':
      argv++;
      if (*argv == NULL) {
        fprintf(ferr, "cylinder: error: -p needs an argument\n");
        exit(1);
      }
      period = strtol(*argv, &end, 10);
      if (*end != '\0' || period <= 0) {
        fprintf(ferr, "cylinder: error: '%s' is not a positive integer\n",
                *argv);
        exit(1);
      }
      PeriodFlag = 1;
      break;
    case 'i':
      Image = 1;
      break;
    case 'e':
      argv++;
      if (*argv == NULL) {
        fprintf(ferr, "cylinder: error: -e needs an argument\n");
        exit(1);
      }
      tend = strtod(*argv, &end);
      if (*end != '\0') {
        fprintf(ferr, "cylinder: error: '%s' is not a number\n", *argv);
        exit(1);
      }
      TendFlag = 1;
      break;
    default:
      fprintf(ferr, "cylinder: error: unknown option '%s'\n", *argv);
      exit(1);
    }
  if (!ReynoldsFlag) {
    fprintf(ferr, "cylinder: error: -r is not set\n");
    exit(1);
  }
  if (!LevelFlag) {
    fprintf(ferr, "cylinder: error: -l is not set\n");
    exit(1);
  }
  if (!PeriodFlag) {
    fprintf(ferr, "cylinder: error: -p is not set\n");
    exit(1);
  }
  if (!TendFlag) {
    fprintf(ferr, "cylinder: error: -e must be set\n");
    exit(1);
  }
  L0 = 4;
  origin((struct _origin){-0.5, -L0 / 2.});
  init_grid(1 << minlevel);
  mu = muv;
  run();
free_solver();}
static int properties_0_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i++);*ip=i;*tp=t;return ret;}      static int properties_0(const int i,const double t,Event *_ev){tracing("properties_0","cylinder.c",0); { foreach_face_stencil(){_stencil_is_face_x(){ {_stencil_val_a(muv.x,0,0,0); _stencil_val(fm.x,0,0,0);     }}end__stencil_is_face_x()_stencil_is_face_y(){ {_stencil_val_a(muv.y,0,0,0); _stencil_val(fm.y,0,0,0);     }}end__stencil_is_face_y()}end_foreach_face_stencil(); if(!is_constant(fm.x)){{foreach_face_generic(){is_face_x(){ val(muv.x,0,0,0) = val(fm.x,0,0,0) * diameter / reynolds;}end_is_face_x()is_face_y(){ val(muv.y,0,0,0) = val(fm.y,0,0,0) * diameter / reynolds;}end_is_face_y()}end_foreach_face_generic();}}else {struct{double x,y;}_const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm); {foreach_face_generic(){is_face_x(){ val(muv.x,0,0,0) = _const_fm.x * diameter / reynolds;}end_is_face_x()is_face_y(){ val(muv.y,0,0,0) = _const_fm.y * diameter / reynolds;}end_is_face_y()}end_foreach_face_generic();}} }{end_tracing("properties_0","cylinder.c",0);return 0;}end_tracing("properties_0","cylinder.c",0);}

static int init_0_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(t = 0);*ip=i;*tp=t;return ret;}      static int init_0(const int i,const double t,Event *_ev){tracing("init_0","cylinder.c",0); {
  scalar  phi=new_vertex_scalar("phi");
  foreach_vertex_stencil() { 
     
       
         
             
    _stencil_val_a(phi,0,0,0);  
  }end_foreach_vertex_stencil();
  {
#line 137
foreach_vertex() {
    double p0;
    p0 = 0.5 - y;
    p0 = min(p0, 0.5 + y);
    p0 = min(p0, sq(x) + sq(y) - sq(diameter / 2));
    val(phi,0,0,0) = p0;
  }end_foreach_vertex();}
  fractions((struct Fractions){phi, cs, fs});
  foreach_stencil () {
    _stencil_val_a(u.x,0,0,0);  
    _stencil_val_a(u.y,0,0,0);  
  }end_foreach_stencil();
  {
#line 145
foreach () {
    val(u.x,0,0,0) = 0;
    val(u.y,0,0,0) = 0;
  }end_foreach();}delete((scalar*)((scalar[]){phi,{-1}}));
}{end_tracing("init_0","cylinder.c",0);return 0;}end_tracing("init_0","cylinder.c",0);}

static int dump_0_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=( t <= tend);*ip=i;*tp=t;return ret;}static int dump_0_expr1(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i++);*ip=i;*tp=t;return ret;}      static int dump_0(const int i,const double t,Event *_ev){tracing("dump_0","cylinder.c",0); {
  static long iframe = 0;
  char png[FILENAME_MAX], htg[FILENAME_MAX];
  scalar  omega=new_scalar("omega"),  m=new_scalar("m");
  FILE *fp;
  long nx, ny;
  coord n, b;

  if (iframe % period == 0) {
    if (pid() == 0)
      fprintf(ferr, "cylinder: %d: %09d %.16e\n", npe(), i, t);
    sprintf(htg, "h.%09ld.htg", iframe);
    vorticity(u, omega);
    output_htg(((scalar[]){p, omega,{-1}}), ((vector[]){u,{{-1},{-1}}}), htg);
    if (Image) {
      foreach_stencil ()
        {_stencil_val_a(m,0,0,0); _stencil_val(cs,0,0,0);   }end_foreach_stencil();
      {
#line 166
foreach ()
        val(m,0,0,0) = val(cs,0,0,0) - 0.5;end_foreach();}
      sprintf(png, "%09ld.ppm", iframe);
      output_ppm((struct OutputPPM){omega, .file = png, .box = {{-0.5, -0.5}, {L0 - 0.5, 0.5}},
                 .min = -5 / diameter, .max = 5 / diameter, .linear = false,
                 .mask = m});
    }
  }
  iframe++;delete((scalar*)((scalar[]){m,omega,{-1}}));
}{end_tracing("dump_0","cylinder.c",0);return 0;}end_tracing("dump_0","cylinder.c",0);}
static int adapt_0_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i++);*ip=i;*tp=t;return ret;}      static int adapt_0(const int i,const double t,Event *_ev){tracing("adapt_0","cylinder.c",0); {
  adapt_wavelet((struct Adapt){((scalar[]){cs, u.x,u.y,{-1}}), (double[]){1e-2, 3e-3, 3e-3}, .maxlevel = level,
                .minlevel = minlevel});
}{end_tracing("adapt_0","cylinder.c",0);return 0;}end_tracing("adapt_0","cylinder.c",0);}
#line 2 "ast/init_solver.h"

static void _init_solver (void)
{
  void init_solver();
  datasize=13*sizeof(double);init_solver();
  embed=new_bid();quadtree_methods();{  event_register((Event){0,1,metric,{metric_expr0},((int *)0),((double *)0),"/u/basilisk/src/embed.h",0,"metric"});
  event_register((Event){0,1,defaults,{defaults_expr0},((int *)0),((double *)0),"/u/basilisk/src/embed.h",0,"defaults"});
  event_register((Event){0,1,defaults_0,{defaults_0_expr0},((int *)0),((double *)0),"/u/basilisk/src/run.h",0,"defaults"});
  event_register((Event){0,1,defaults_1,{defaults_1_expr0},((int *)0),((double *)0),"/u/basilisk/src/navier-stokes/centered.h",0,"defaults"});
  event_register((Event){0,1,default_display,{default_display_expr0},((int *)0),((double *)0),"/u/basilisk/src/navier-stokes/centered.h",0,"default_display"});
  event_register((Event){0,1,init,{init_expr0},((int *)0),((double *)0),"/u/basilisk/src/navier-stokes/centered.h",0,"init"});
  event_register((Event){0,1,init_0,{init_0_expr0},((int *)0),((double *)0),"cylinder.c",0,"init"});
  event_register((Event){0,2,dump_0,{dump_0_expr0,dump_0_expr1},((int *)0),((double *)0),"cylinder.c",0,"dump"});


    

    init_const_vector((vector){{_NVARMAX+0},{_NVARMAX+1}},"zerof",(double[]) {0.,0.,0.});
  init_const_vector((vector){{_NVARMAX+2},{_NVARMAX+3}},"unityf",(double[]) {1.,1.,1.});
  init_const_scalar((scalar){_NVARMAX+4},"unity", 1.);
  init_const_scalar((scalar){_NVARMAX+5},"zeroc", 0.);
  init_scalar((scalar){0},"cs");
  init_face_vector((vector){{1},{2}},"fs");
  event_register((Event){0,1,cleanup,{cleanup_expr0},((int *)0),((double *)0),"/u/basilisk/src/run.h",0,"cleanup"});
  init_scalar((scalar){3},"p");
  init_vector((vector){{4},{5}},"u");
  init_vector((vector){{6},{7}},"g");
  init_scalar((scalar){8},"pf");
  init_face_vector((vector){{9},{10}},"uf");
  event_register((Event){0,1,set_dtmax,{set_dtmax_expr0},((int *)0),((double *)0),"/u/basilisk/src/navier-stokes/centered.h",0,"set_dtmax"});
  event_register((Event){0,1,stability,{stability_expr0},((int *)0),((double *)0),"/u/basilisk/src/navier-stokes/centered.h",0,"stability"});
  event_register((Event){0,1,vof,{vof_expr0},((int *)0),((double *)0),"/u/basilisk/src/navier-stokes/centered.h",0,"vof"});
  event_register((Event){0,1,tracer_advection,{tracer_advection_expr0},((int *)0),((double *)0),"/u/basilisk/src/navier-stokes/centered.h",0,"tracer_advection"});
  event_register((Event){0,1,tracer_diffusion,{tracer_diffusion_expr0},((int *)0),((double *)0),"/u/basilisk/src/navier-stokes/centered.h",0,"tracer_diffusion"});
  event_register((Event){0,1,properties,{properties_expr0},((int *)0),((double *)0),"/u/basilisk/src/navier-stokes/centered.h",0,"properties"});
  event_register((Event){0,1,advection_term,{advection_term_expr0},((int *)0),((double *)0),"/u/basilisk/src/navier-stokes/centered.h",0,"advection_term"});
  event_register((Event){0,1,viscous_term,{viscous_term_expr0},((int *)0),((double *)0),"/u/basilisk/src/navier-stokes/centered.h",0,"viscous_term"});
  event_register((Event){0,1,acceleration,{acceleration_expr0},((int *)0),((double *)0),"/u/basilisk/src/navier-stokes/centered.h",0,"acceleration"});
  event_register((Event){0,1,projection,{projection_expr0},((int *)0),((double *)0),"/u/basilisk/src/navier-stokes/centered.h",0,"projection"});
  event_register((Event){0,1,end_timestep,{end_timestep_expr0},((int *)0),((double *)0),"/u/basilisk/src/navier-stokes/centered.h",0,"end_timestep"});
  event_register((Event){0,1,adapt,{adapt_expr0},((int *)0),((double *)0),"/u/basilisk/src/navier-stokes/centered.h",0,"adapt"});
  init_face_vector((vector){{11},{12}},"muv");
  event_register((Event){0,1,properties_0,{properties_0_expr0},((int *)0),((double *)0),"cylinder.c",0,"properties"});
  event_register((Event){0,1,adapt_0,{adapt_0_expr0},((int *)0),((double *)0),"cylinder.c",0,"adapt"});

#line 13
}  _attribute[p.i].dirty=1,_attribute[p.i].boundary[right]=_boundary0,_attribute[p.i].boundary_homogeneous[right]=_boundary0_homogeneous;
  _attribute[p.i].dirty=1,_attribute[p.i].boundary[left]=_boundary1,_attribute[p.i].boundary_homogeneous[left]=_boundary1_homogeneous;
  _attribute[p.i].dirty=1,_attribute[p.i].boundary[top]=_boundary2,_attribute[p.i].boundary_homogeneous[top]=_boundary2_homogeneous;
  _attribute[p.i].dirty=1,_attribute[p.i].boundary[bottom]=_boundary3,_attribute[p.i].boundary_homogeneous[bottom]=_boundary3_homogeneous;
  _attribute[u.x.i].dirty=1,_attribute[u.x.i].boundary[left]=_boundary4,_attribute[u.x.i].boundary_homogeneous[left]=_boundary4_homogeneous;
  _attribute[p.i].dirty=1,_attribute[p.i].boundary[left]=_boundary5,_attribute[p.i].boundary_homogeneous[left]=_boundary5_homogeneous;
  _attribute[pf.i].dirty=1,_attribute[pf.i].boundary[left]=_boundary6,_attribute[pf.i].boundary_homogeneous[left]=_boundary6_homogeneous;
  _attribute[u.x.i].dirty=1,_attribute[u.x.i].boundary[right]=_boundary7,_attribute[u.x.i].boundary_homogeneous[right]=_boundary7_homogeneous;
  _attribute[p.i].dirty=1,_attribute[p.i].boundary[right]=_boundary8,_attribute[p.i].boundary_homogeneous[right]=_boundary8_homogeneous;
  _attribute[pf.i].dirty=1,_attribute[pf.i].boundary[right]=_boundary9,_attribute[pf.i].boundary_homogeneous[right]=_boundary9_homogeneous;
  _attribute[u.x.i].dirty=1,_attribute[u.x.i].boundary[embed]=_boundary10,_attribute[u.x.i].boundary_homogeneous[embed]=_boundary10_homogeneous;
  _attribute[u.y.i].dirty=1,_attribute[u.y.i].boundary[embed]=_boundary11,_attribute[u.y.i].boundary_homogeneous[embed]=_boundary11_homogeneous;

  
#line 14
set_fpe();
}
