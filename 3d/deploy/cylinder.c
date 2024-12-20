#line 0 "cylinder-cpp.c"
#line 0 "<built-in>"
#line 0 "<command-line>"
#line 1 "/usr/include/stdc-predef.h"
#line 0 "<command-line>"
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



#line 1 "/home/lisergey/basilisk/src/common.h"
#line 1 "/home/lisergey/basilisk/src/ast/std/stdlib.h"
#include <stdlib.h>
#line 2 "/home/lisergey/basilisk/src/common.h"
#line 1 "/home/lisergey/basilisk/src/ast/std/stdio.h"
#include <stdio.h>
#line 3 "/home/lisergey/basilisk/src/common.h"
#line 1 "/home/lisergey/basilisk/src/ast/std/stddef.h"
#include <stddef.h>
#line 4 "/home/lisergey/basilisk/src/common.h"
#line 1 "/home/lisergey/basilisk/src/ast/std/stdbool.h"
#include <stdbool.h>
#line 5 "/home/lisergey/basilisk/src/common.h"
#line 1 "/home/lisergey/basilisk/src/ast/std/stdarg.h"
#include <stdarg.h>
#line 6 "/home/lisergey/basilisk/src/common.h"
#line 1 "/home/lisergey/basilisk/src/ast/std/string.h"
#include <string.h>
#line 7 "/home/lisergey/basilisk/src/common.h"
#line 1 "/home/lisergey/basilisk/src/ast/std/float.h"
#include <float.h>
#line 8 "/home/lisergey/basilisk/src/common.h"
#line 1 "/home/lisergey/basilisk/src/ast/std/limits.h"
#include <limits.h>
#line 9 "/home/lisergey/basilisk/src/common.h"
#line 1 "/home/lisergey/basilisk/src/ast/std/math.h"
#include <math.h>
#line 10 "/home/lisergey/basilisk/src/common.h"
#line 1 "/home/lisergey/basilisk/src/ast/std/time.h"
#include <time.h>
#line 11 "/home/lisergey/basilisk/src/common.h"
#line 1 "/home/lisergey/basilisk/src/ast/std/sys/time.h"
#include <sys/time.h>
#line 12 "/home/lisergey/basilisk/src/common.h"
#line 1 "/home/lisergey/basilisk/src/ast/std/sys/resource.h"
#include <sys/resource.h>
#line 13 "/home/lisergey/basilisk/src/common.h"

#if _OPENMP
# include <omp.h>
# define OMP(x) _Pragma(#x)
#elif 1

# define OMP(x)

# include <mpi.h>
static int mpi_rank, mpi_npe;
# define tid() mpi_rank
# define pid() mpi_rank
# define npe() mpi_npe

#else

# define OMP(x)

#endif
#line 49 "/home/lisergey/basilisk/src/common.h"
#define _NVARMAX 65536
#define is_constant(v) ((v).i >= _NVARMAX)
#define constant(v) (is_constant(v) ? _constant[(v).i - _NVARMAX] : 1e30)

#define max(a,b) ((a) > (b) ? (a) : (b))
#define min(a,b) ((a) < (b) ? (a) : (b))
#define sq(x) ((x)*(x))
#define cube(x) ((x)*(x)*(x))
#define sign(x) ((x) > 0 ? 1 : -1)
#define sign2(x) ((x) > 0 ? 1 : (x) < 0 ? -1 : 0)
#define noise() (1. - 2.*rand()/(double)RAND_MAX)
#define clamp(x,a,b) ((x) < (a) ? (a) : (x) > (b) ? (b) : (x))

#define unmap(x,y)

#define trash(x)


#define systderr stderr
#define systdout stdout

#if 1
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

#line 81

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
#line 105 "/home/lisergey/basilisk/src/common.h"
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
  if (!(d != NULL)) qassert ("/home/lisergey/basilisk/src/common.h", 180, "d != NULL");
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
  Array * a = ((Array *) pmalloc ((1)*sizeof(Array),__func__,__FILE__,__LINE__));
  a->p = NULL;
  a->max = a->len = 0;
  return a;
}

void array_free (Array * a)
{
  pfree (a->p,__func__,__FILE__,__LINE__);
  pfree (a,__func__,__FILE__,__LINE__);
}

void array_append (Array * a, void * elem, size_t size)
{
  if (a->len + size >= a->max) {
    a->max += max (size, 4096);
    a->p = prealloc (a->p, a->max,__func__,__FILE__,__LINE__);
  }
  memcpy (((char *)a->p) + a->len, elem, size);
  a->len += size;
}

void * array_shrink (Array * a)
{
  void * p = prealloc (a->p, a->len,__func__,__FILE__,__LINE__);
  pfree (a,__func__,__FILE__,__LINE__);
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
  char * s = pstrdup (func,__func__,__FILE__,__LINE__);
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
  if (!(t->stack.len > 0)) qassert ("/home/lisergey/basilisk/src/common.h", 452, "t->stack.len > 0");
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
    pfree (*func,__func__,__FILE__,__LINE__);
  pfree (t->index.p,__func__,__FILE__,__LINE__);
  pfree (t->stack.p,__func__,__FILE__,__LINE__);
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
#if 1
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
    TraceIndex t = {pstrdup(func,__func__,__FILE__,__LINE__), pstrdup(file,__func__,__FILE__,__LINE__), line, 1, total, self};
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
  if (!(Trace.stack.len >= 2*sizeof(double))) qassert ("/home/lisergey/basilisk/src/common.h", 556, "Trace.stack.len >= 2*sizeof(double)");
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

#if 1
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
#if 1
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
#if 1
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
    pfree (t->func,__func__,__FILE__,__LINE__), pfree (t->file,__func__,__FILE__,__LINE__);

  pfree (Trace.index.p,__func__,__FILE__,__LINE__);
  Trace.index.p = NULL;
  Trace.index.len = Trace.index.max = 0;

  pfree (Trace.stack.p,__func__,__FILE__,__LINE__);
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

#elif 1

static bool in_prof = false;
static double prof_start, _prof;
#define prof_start(name)\
  if (!(!in_prof)) qassert ("/home/lisergey/basilisk/src/common.h", 666, "!in_prof"); in_prof = true;\
  prof_start = MPI_Wtime();\

#line 668

#define prof_stop()\
  if (!(in_prof)) qassert ("/home/lisergey/basilisk/src/common.h", 670, "in_prof"); in_prof = false;\
  _prof = MPI_Wtime();\
  mpi_time += _prof - prof_start;\

#line 673


#if FAKE_MPI
#define mpi_all_reduce(v,type,op)
#define mpi_all_reduce_array(v,type,op,elem)
#else
     
int mpi_all_reduce0 (void *sendbuf, void *recvbuf, int count,
       MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
{tracing("mpi_all_reduce0","/home/lisergey/basilisk/src/common.h",680);
  { int _ret= MPI_Allreduce (sendbuf, recvbuf, count, datatype, op, comm);end_tracing("mpi_all_reduce0","/home/lisergey/basilisk/src/common.h",683);return _ret;}
end_tracing("mpi_all_reduce0","/home/lisergey/basilisk/src/common.h",684);}
#define mpi_all_reduce(v,type,op) {\
  prof_start ("mpi_all_reduce");\
  union { int a; float b; double c;} global;\
  mpi_all_reduce0 (&(v), &global, 1, type, op, MPI_COMM_WORLD);\
  memcpy (&(v), &global, sizeof (v));\
  prof_stop();\
}\

#line 692

#define mpi_all_reduce_array(v,type,op,elem) {\
  prof_start ("mpi_all_reduce");\
  type * global = malloc ((elem)*sizeof(type)), * tmp = malloc ((elem)*sizeof(type));\
  for (int i = 0; i < elem; i++)\
    tmp[i] = (v)[i];\
  MPI_Datatype datatype;\
  if (!strcmp(#type, "double")) datatype = MPI_DOUBLE;\
  else if (!strcmp(#type, "int")) datatype = MPI_INT;\
  else if (!strcmp(#type, "long")) datatype = MPI_LONG;\
  else if (!strcmp(#type, "bool")) datatype = MPI_C_BOOL;\
  else if (!strcmp(#type, "unsigned char")) datatype = MPI_UNSIGNED_CHAR;\
  else {\
    fprintf (stderr, "unknown reduction type '%s'\n", #type);\
    fflush (stderr);\
    abort();\
  }\
  mpi_all_reduce0 (tmp, global, elem, datatype, op, MPI_COMM_WORLD);\
  for (int i = 0; i < elem; i++)\
    (v)[i] = global[i];\
  free (global), free (tmp);\
  prof_stop();\
}\

#line 715


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
#line 826 "/home/lisergey/basilisk/src/common.h"
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
  if (!(sizeof (int64_t) == sizeof (double))) qassert ("/home/lisergey/basilisk/src/common.h", 836, "sizeof (int64_t) == sizeof (double)");
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




int N = 16;


typedef struct { int i; } scalar;

typedef struct {
  scalar x;

  scalar y;


  scalar z;

} vector;

typedef struct {
  scalar * x;

  scalar * y;


  scalar * z;

} vectorl;

typedef struct {
  vector x;

  vector y;


  vector z;

} tensor;

struct { int x, y, z; } Period = {false, false, false};

typedef struct {
  double x, y, z;
} coord;

OMP(omp declare reduction (+ : coord :
      omp_out.x += omp_in.x,
      omp_out.y += omp_in.y,
      omp_out.z += omp_in.z))
#line 920 "/home/lisergey/basilisk/src/common.h"
void normalize (coord * n)
{
  double norm = 0.;
  
    norm += sq(n->x);
    
#line 924
norm += sq(n->y);
    
#line 924
norm += sq(n->z);
  norm = sqrt(norm);
  
    n->x /= norm;
    
#line 927
n->y /= norm;
    
#line 927
n->z /= norm;
}

void origin (double x, double y, double z) {
  X0 = x; Y0 = y; Z0 = z;
}

void size (double L) {
  L0 = L;
}

double zero (double s0, double s1, double s2) { return 0.; }
#line 947 "/home/lisergey/basilisk/src/common.h"
  enum { right, left, top, bottom, front, back };

int nboundary = 2*3;



#define _dirichlet(expr, ...) (2.*(expr) - val(_s,0,0,0))
#define _dirichlet_homogeneous(...) (- val(_s,0,0,0))
#define _dirichlet_face(expr,...) (expr)
#define _dirichlet_face_homogeneous(...) (0.)
#define _neumann(expr,...) (Delta*(expr) + val(_s,0,0,0))
#define _neumann_homogeneous(...) (val(_s,0,0,0))

double * _constant = NULL;
size_t datasize = 0;
typedef struct _Point Point;

#line 1 "/home/lisergey/basilisk/src/grid/boundaries.h"


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
  boundaries = (Boundary * *) prealloc (boundaries, (len + 2)*sizeof(Boundary *),__func__,__FILE__,__LINE__);
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
      pfree (b,__func__,__FILE__,__LINE__);
  pfree (boundaries,__func__,__FILE__,__LINE__);
  boundaries = NULL;
}
#line 47 "/home/lisergey/basilisk/src/grid/boundaries.h"
typedef struct {
  Boundary parent;
  int d;
} BoxBoundary;
#line 965 "/home/lisergey/basilisk/src/common.h"



typedef struct {
  double (** boundary) (Point, Point, scalar, void *);
  double (** boundary_homogeneous) (Point, Point, scalar, void *);
  double (* gradient) (double, double, double);
  void (* delete) (scalar);
  char * name;
  struct {
    int x;

    int y;


    int z;

  } d;
  vector v;
  int face;
  bool nodump, freed;
  int block;
  scalar * depends;

  
#line 19 "/home/lisergey/basilisk/src/grid/stencils.h"
bool input, output;
  int width;
  int dirty;
  
#line 18 "/home/lisergey/basilisk/src/grid/multigrid-common.h"
void (* prolongation) (Point, scalar);
  void (* restriction) (Point, scalar);
  
#line 9 "/home/lisergey/basilisk/src/grid/tree-common.h"
void (* refine) (Point, scalar);
  
#line 97
void (* coarsen) (Point, scalar);
  
#line 81 "/home/lisergey/basilisk/src/fractions.h"
vector n;

#line 988 "/home/lisergey/basilisk/src/common.h"
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
  list = (scalar *) prealloc (list, (len + 2)*sizeof(scalar),__func__,__FILE__,__LINE__);
  list[len] = s;
  list[len + 1].i = -1;
  return list;
}

scalar * list_prepend (scalar * list, scalar s)
{
  int len = list_len (list);
  list = (scalar *) prealloc (list, (len + 2)*sizeof(scalar),__func__,__FILE__,__LINE__);
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
  list = (vector *) prealloc (list, (len + 2)*sizeof(vector),__func__,__FILE__,__LINE__);
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
      
#line 1089
if (w.y.i != v.y.i)
 id = false;
      
#line 1089
if (w.z.i != v.z.i)
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
      if (!(s->i >= 0)) qassert ("/home/lisergey/basilisk/src/common.h", 1112, "s->i >= 0");
      v.x = *s++;
    } 
#line 1111
{
      if (!(s->i >= 0)) qassert ("/home/lisergey/basilisk/src/common.h", 1112, "s->i >= 0");
      v.y = *s++;
    } 
#line 1111
{
      if (!(s->i >= 0)) qassert ("/home/lisergey/basilisk/src/common.h", 1112, "s->i >= 0");
      v.z = *s++;
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
  list = (tensor *) prealloc (list, (len + 2)*sizeof(tensor),__func__,__FILE__,__LINE__);
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
      if (!(v->x.i >= 0)) qassert ("/home/lisergey/basilisk/src/common.h", 1143, "v->x.i >= 0");
      t.x = *v++;
    } 
#line 1142
{
      if (!(v->y.i >= 0)) qassert ("/home/lisergey/basilisk/src/common.h", 1143, "v->x.i >= 0");
      t.y = *v++;
    } 
#line 1142
{
      if (!(v->z.i >= 0)) qassert ("/home/lisergey/basilisk/src/common.h", 1143, "v->x.i >= 0");
      t.z = *v++;
    }
    list = tensors_append (list, t);
  }
  return list;
}

static inline bool is_vertex_scalar (scalar s)
{
  
    if (_attribute[s.i].d.x != -1)
      return false;
    
#line 1154
if (_attribute[s.i].d.y != -1)
      return false;
    
#line 1154
if (_attribute[s.i].d.z != -1)
      return false;
  return true;
}

scalar * all = NULL;
scalar * baseblock = NULL;



scalar (* init_scalar) (scalar, const char *);
scalar (* init_vertex_scalar) (scalar, const char *);
vector (* init_vector) (vector, const char *);
vector (* init_face_vector) (vector, const char *);
tensor (* init_tensor) (tensor, const char *);
void (* scalar_clone) (scalar, scalar);





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
  Events = pmalloc (sizeof (Event),__func__,__FILE__,__LINE__);
  Events[0].last = 1;
  _attribute = pcalloc (datasize/sizeof(double), sizeof (_Attributes),__func__,__FILE__,__LINE__);
  int n = datasize/sizeof(double);
  all = (scalar *) pmalloc (sizeof (scalar)*(n + 1),__func__,__FILE__,__LINE__);
  baseblock = (scalar *) pmalloc (sizeof (scalar)*(n + 1),__func__,__FILE__,__LINE__);
  for (int i = 0; i < n; i++)
    baseblock[i].i = all[i].i = i;
  baseblock[n].i = all[n].i = -1;
#if _CADNA
  cadna_init (-1);
#endif
#if 1
  mpi_init();
#elif MTRACE == 1
  char * etrace = getenv ("MTRACE");
  pmtrace.fp = fopen (etrace ? etrace : "mtrace", "w");
  pmtrace.fname = systrdup (etrace ? etrace : "mtrace");
#endif
}



#if 1
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
#if 1
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



const vector zerof = {{_NVARMAX+0},{_NVARMAX+1},{_NVARMAX+2}};
const vector unityf = {{_NVARMAX+3},{_NVARMAX+4},{_NVARMAX+5}};
const scalar unity = {_NVARMAX+6};
const scalar zeroc = {_NVARMAX+7};



const vector unityf0 = {{_NVARMAX+8},{_NVARMAX+9},{_NVARMAX+10}};
const scalar unity0 = {_NVARMAX+11};
        vector fm = {{_NVARMAX+8},{_NVARMAX+9},{_NVARMAX+10}};
        scalar cm = {_NVARMAX+11};
#line 1281 "/home/lisergey/basilisk/src/common.h"
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
    qpopen_pipes = (FILE * *) prealloc (qpopen_pipes, (n + 2)*sizeof(FILE *),__func__,__FILE__,__LINE__);
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
  pfree (qpopen_pipes,__func__,__FILE__,__LINE__);
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
  void ** m = ((void * *) pmalloc ((n)*sizeof(void *),__func__,__FILE__,__LINE__));
  char * a = ((char *) pmalloc ((n*p*size)*sizeof(char),__func__,__FILE__,__LINE__));
  for (int i = 0; i < n; i++)
    m[i] = a + i*p*size;
  return m;
}

double matrix_inverse (double ** m, int n, double pivmin)
{
  int indxc[n], indxr[n], ipiv[n];
  int i, icol = 0, irow = 0, j, k, l, ll;
  double big, dum, pivinv, minpiv = 1e30;

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
  pfree (((void **) m)[0],__func__,__FILE__,__LINE__);
  pfree (m,__func__,__FILE__,__LINE__);
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

static void free_display_defaults() {
  pfree (display_defaults,__func__,__FILE__,__LINE__);
}

void display (const char * commands, bool overwrite)
{
  if (display_defaults == NULL)
    free_solver_func_add (free_display_defaults);
  if (overwrite) {
    pfree (display_defaults,__func__,__FILE__,__LINE__);
    display_defaults = pmalloc (strlen(commands) + 2,__func__,__FILE__,__LINE__);
    strcpy (display_defaults, "@");
    strcat (display_defaults, commands);
  }
  else {
    if (!display_defaults)
      display_defaults = pstrdup ("@",__func__,__FILE__,__LINE__);
    display_defaults =
      prealloc (display_defaults,
        strlen(display_defaults) + strlen(commands) + 1,__func__,__FILE__,__LINE__);
    strcat (display_defaults, commands);
  }
}



typedef struct {
  double x;

  double y;


  double z;

} _coord;





#line 1 "/home/lisergey/basilisk/src/grid/stencils.h"
#line 17 "/home/lisergey/basilisk/src/grid/stencils.h"










typedef struct {
  char * name;
  char * type;
  void * pointer;
  int * dimensions;
  int nd;
  char reduct;
  scalar data;
} NonLocal;

typedef struct {
  const char * fname;
  int line;
  int first;
  int face;
  bool vertex;
  int parallel;
  scalar * listc;
  vectorl listf;
  scalar * dirty;
  void * data;
} ForeachData;


#define foreach_stencil(...) {\
  static int _first = 1.;\
  ForeachData _loop = {\
    .fname = __FILE__, .line = __LINE__, .first = _first\
  };\
  if (baseblock) for (scalar s = baseblock[0], * i = baseblock;\
  s.i >= 0; i++, s = *i) {\
    _attribute[s.i].input = _attribute[s.i].output = false;\
    _attribute[s.i].width = 0;\
  }\
  int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);\
  Point point = {0}; NOT_UNUSED (point);\

#line 63


#define end_foreach_stencil()\
  check_stencil (&_loop);\
  boundary_stencil (&_loop);\
  _first = 0;\
}\

#line 70


#define foreach_vertex_stencil(...) foreach_stencil(__VA_ARGS__) _loop.vertex = true;
#define end_foreach_vertex_stencil() end_foreach_stencil()

#define foreach_face_stencil(...) foreach_stencil(__VA_ARGS__)
#define end_foreach_face_stencil() end_foreach_stencil()

#define foreach_visible_stencil(...) foreach_stencil(__VA_ARGS__)
#define end_foreach_visible_stencil() end_foreach_stencil()

#define foreach_point_stencil(...) foreach_stencil(__VA_ARGS__)
#define end_foreach_point_stencil() end_foreach_stencil()

#define foreach_region_stencil(...) foreach_stencil(__VA_ARGS__)
#define end_foreach_region_stencil() end_foreach_stencil()

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
  stencil_val (point, a, _i, _j, _k, __FILE__, __LINE__, false)\

#line 101

#define _stencil_val_o(a,_i,_j,_k)\
  stencil_val (point, a, _i, _j, _k, __FILE__, __LINE__, true)\

#line 104

#define _stencil_val_a(a,_i,_j,_k)\
  stencil_val_a (point, a, _i, _j, _k, false, __FILE__, __LINE__)\

#line 107

#define _stencil_val_r(a,_i,_j,_k)\
  stencil_val_a (point, a, _i, _j, _k, true, __FILE__, __LINE__)\

#line 110


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







void check_stencil (ForeachData * loop)
{
  loop->listf = (vectorl){NULL};




  {scalar*_i=(scalar*)( baseblock);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {
    bool write = _attribute[s.i].output, read = _attribute[s.i].input;




    {





      if (read && scalar_is_dirty (s)) {





 if (_attribute[s.i].face) {
   if (_attribute[s.i].width > 0)
     loop->listc = list_append (loop->listc, s);
   else if (!write) {
     scalar sn = _attribute[s.i].v.x.i >= 0 ? _attribute[s.i].v.x : s;
     
       if (_attribute[s.i].v.x.i == s.i) {




  if (_attribute[sn.i].boundary[left] || _attribute[sn.i].boundary[right])
    loop->listc = list_append (loop->listc, s);
  else if (_attribute[s.i].dirty != 2)
    loop->listf.x = list_append (loop->listf.x, s);
       }
       
#line 214
if (_attribute[s.i].v.y.i == s.i) {




  if (_attribute[sn.i].boundary[bottom] || _attribute[sn.i].boundary[top])
    loop->listc = list_append (loop->listc, s);
  else if (_attribute[s.i].dirty != 2)
    loop->listf.y = list_append (loop->listf.y, s);
       }
       
#line 214
if (_attribute[s.i].v.z.i == s.i) {




  if (_attribute[sn.i].boundary[back] || _attribute[sn.i].boundary[front])
    loop->listc = list_append (loop->listc, s);
  else if (_attribute[s.i].dirty != 2)
    loop->listf.z = list_append (loop->listf.z, s);
       }
   }
 }





 else if (_attribute[s.i].width > 0)
   loop->listc = list_append (loop->listc, s);
      }





      if (write) {
 if (3 > 1 && !loop->vertex && loop->first) {
   bool vertex = true;
   
     if (_attribute[s.i].d.x != -1)
       vertex = false;
     
#line 243
if (_attribute[s.i].d.y != -1)
       vertex = false;
     
#line 243
if (_attribute[s.i].d.z != -1)
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
#line 261
{
       if (loop->face == d) {
  _attribute[s.i].face = 2, _attribute[s.i].v.y.i = s.i;
  _attribute[s.i].boundary[bottom] = _attribute[s.i].boundary[top] = NULL;





       }
       d *= 2, i++;
     } 
#line 261
{
       if (loop->face == d) {
  _attribute[s.i].face = 2, _attribute[s.i].v.z.i = s.i;
  _attribute[s.i].boundary[back] = _attribute[s.i].boundary[front] = NULL;





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
       name = pstrdup (_attribute[s.i].name,__func__,__FILE__,__LINE__);
       char * s = name + strlen(name) - 1;
       while (s != name && *s != '.') s--;
       if (s != name) *s = '\0';
     }
     struct { int x, y, z; } input, output;
     vector v = _attribute[s.i].v;

     
       input.x = _attribute[v.x.i].input, output.x = _attribute[v.x.i].output;
       
#line 291
input.y = _attribute[v.y.i].input, output.y = _attribute[v.y.i].output;
       
#line 291
input.z = _attribute[v.z.i].input, output.z = _attribute[v.z.i].output;

     init_face_vector (v, name);


     
       _attribute[v.x.i].input = input.x, _attribute[v.x.i].output = output.x;
       
#line 297
_attribute[v.y.i].input = input.y, _attribute[v.y.i].output = output.y;
       
#line 297
_attribute[v.z.i].input = input.z, _attribute[v.z.i].output = output.z;





     pfree (name,__func__,__FILE__,__LINE__);
   }
 }
 else if (loop->vertex) {
   bool vertex = true;
   
     if (_attribute[s.i].d.x != -1)
       vertex = false;
     
#line 309
if (_attribute[s.i].d.y != -1)
       vertex = false;
     
#line 309
if (_attribute[s.i].d.z != -1)
       vertex = false;
   if (!vertex) {
     char * name = NULL;
     if (_attribute[s.i].name) name = pstrdup (_attribute[s.i].name,__func__,__FILE__,__LINE__);
     init_vertex_scalar (s, name);
     
       _attribute[s.i].v.x.i = -1;
       
#line 316
_attribute[s.i].v.y.i = -1;
       
#line 316
_attribute[s.i].v.z.i = -1;




     pfree (name,__func__,__FILE__,__LINE__);
   }
 }





 loop->dirty = list_append (loop->dirty, s);
 {scalar*_i=(scalar*)( baseblock);if(_i)for(scalar d=*_i;(&d)->i>=0;d=*++_i){
   if (scalar_depends_from (d, s))
     loop->dirty = list_append (loop->dirty, d);}}
      }
    }
  }}}
}




void boundary_stencil (ForeachData * loop)
{
  bool flux = false;
  
    if (loop->listf.x)
      flux = true;
    
#line 345
if (loop->listf.y)
      flux = true;
    
#line 345
if (loop->listf.z)
      flux = true;
  if (flux) {
#line 360 "/home/lisergey/basilisk/src/grid/stencils.h"
    boundary_face (loop->listf);
    
      pfree (loop->listf.x,__func__,__FILE__,__LINE__), loop->listf.x = NULL;
      
#line 362
pfree (loop->listf.y,__func__,__FILE__,__LINE__), loop->listf.y = NULL;
      
#line 362
pfree (loop->listf.z,__func__,__FILE__,__LINE__), loop->listf.z = NULL;
  }




  if (loop->listc) {






    boundary_internal (loop->listc, loop->fname, loop->line);
    pfree (loop->listc,__func__,__FILE__,__LINE__), loop->listc = NULL;
  }





  if (loop->dirty) {






    {scalar*_i=(scalar*)( loop->dirty);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      _attribute[s.i].dirty = true;}}
    pfree (loop->dirty,__func__,__FILE__,__LINE__), loop->dirty = NULL;
  }
}
#line 1462 "/home/lisergey/basilisk/src/common.h"





typedef struct {
  coord x, y, z;
} mat3;

OMP(omp declare reduction (+ : mat3 :
      omp_out.x.x += omp_in.x.x,
      omp_out.x.y += omp_in.x.y,
      omp_out.x.z += omp_in.x.z,
      omp_out.y.x += omp_in.y.x,
      omp_out.y.y += omp_in.y.y,
      omp_out.y.z += omp_in.y.z,
      omp_out.z.x += omp_in.z.x,
      omp_out.z.y += omp_in.z.y,
      omp_out.z.z += omp_in.z.z
      ))
#line 14 "cylinder-cpp.c"
#line 1 "cylinder.c"
#include <float.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#line 1 "grid/octree.h"
#line 1 "/home/lisergey/basilisk/src/grid/octree.h"


#line 1 "grid/tree.h"
#line 1 "/home/lisergey/basilisk/src/grid/tree.h"
#line 1 "grid/mempool.h"
#line 1 "/home/lisergey/basilisk/src/grid/mempool.h"





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

  if (!(poolsize % 8 == 0)) qassert ("/home/lisergey/basilisk/src/grid/mempool.h", 26, "poolsize % 8 == 0");
  if (!(size >= sizeof(FreeBlock))) qassert ("/home/lisergey/basilisk/src/grid/mempool.h", 27, "size >= sizeof(FreeBlock)");


  poolsize = min(1 << 20, poolsize + sizeof(Pool));
  Mempool * m = ((Mempool *) pcalloc (1, sizeof(Mempool),__func__,__FILE__,__LINE__));
  m->poolsize = poolsize;
  m->size = size;
  return m;
}

void mempool_destroy (Mempool * m)
{
  Pool * p = m->pool;
  while (p) {
    Pool * next = p->next;
    pfree (p,__func__,__FILE__,__LINE__);
    p = next;
  }
  pfree (m,__func__,__FILE__,__LINE__);
}

void * mempool_alloc (Mempool * m)
{
  if (!m->first) {

    Pool * p = (Pool *) pmalloc (m->poolsize,__func__,__FILE__,__LINE__);
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
#line 2 "/home/lisergey/basilisk/src/grid/tree.h"




#line 1 "grid/memindex/range.h"
#line 1 "/home/lisergey/basilisk/src/grid/memindex/range.h"
#line 15 "/home/lisergey/basilisk/src/grid/memindex/range.h"
typedef struct {
  void ** p;
  int size;
} Memalloc;

typedef struct {
  int start, end;
} Memrange;
#line 34 "/home/lisergey/basilisk/src/grid/memindex/range.h"
void memrange_alloc (Memrange * r, Memalloc * mem, int i)
{
  if (r->start == r->end) {
    r->start = i;
    r->end = i + 1;
    for (Memalloc * m = mem; m->p; m++) {
      *m->p = pcalloc (1, m->size,__func__,__FILE__,__LINE__);
      *m->p = (char *)(*m->p) - i*m->size;
    }
  }
  else if (i >= r->end) {
    for (Memalloc * m = mem; m->p; m++) {
      *m->p = prealloc ((char *)(*m->p) + r->start*m->size,
         m->size*(i + 1 - r->start),__func__,__FILE__,__LINE__);
      *m->p = (char *)(*m->p) - r->start*m->size;
      memset ((char *)(*m->p) + r->end*m->size, 0, (i - r->end + 1)*m->size);
    }
    r->end = i + 1;
  }
  else if (i < r->start) {
    for (Memalloc * m = mem; m->p; m++) {
      *m->p = prealloc ((char *)(*m->p) + r->start*m->size, m->size*(r->end - i),__func__,__FILE__,__LINE__);
      memmove ((char *)(*m->p) + (r->start - i)*m->size, *m->p,
        m->size*(r->end - r->start));
      memset ((char *)(*m->p), 0, (r->start - i)*m->size);
      *m->p = (char *)(*m->p) - i*m->size;
    }
    r->start = i;
  }
}
#line 73 "/home/lisergey/basilisk/src/grid/memindex/range.h"
bool memrange_free (Memrange * r, Memalloc * mem, int i)
{
  if (i == r->start) {
    if (i == r->end - 1) {
      for (Memalloc * m = mem; m->p; m++) {
 pfree ((char *)(*m->p) + r->start*m->size,__func__,__FILE__,__LINE__);
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
    m->size*(r->end - i),__func__,__FILE__,__LINE__);
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
         m->size*(r->end - r->start),__func__,__FILE__,__LINE__);
      *m->p = (char *)(*m->p) - r->start*m->size;
    }
  }
  else {
    if (!(i > r->start && i < r->end)) qassert ("/home/lisergey/basilisk/src/grid/memindex/range.h", 108, "i > r->start && i < r->end");
    for (Memalloc * m = mem; m->p; m++)
      memset ((char *)(*m->p) + i*m->size, 0, m->size);
  }
  return false;
}







struct _Memindex {
  Memrange r1;

  Memrange * r2;

  Memrange ** r3;







  char **** b;

};
#line 171 "/home/lisergey/basilisk/src/grid/memindex/range.h"
struct _Memindex * mem_new (int len)
{
  struct _Memindex * m = pcalloc (1, sizeof (struct _Memindex),__func__,__FILE__,__LINE__);
  return m;
}





void mem_destroy (struct _Memindex * m, int len)
{

  for (int i = m->r1.start; i < m->r1.end; i++)
    if (m->b[i]) {

      for (int j = m->r2[i].start; j < m->r2[i].end; j++)
 if (m->b[i][j])
   pfree (m->b[i][j] + m->r3[i][j].start,__func__,__FILE__,__LINE__);
      pfree (m->r3[i] + m->r2[i].start,__func__,__FILE__,__LINE__);

      pfree (m->b[i] + m->r2[i].start,__func__,__FILE__,__LINE__);
    }
  if (m->b) {
    pfree (m->r2 + m->r1.start,__func__,__FILE__,__LINE__);

    pfree (m->r3 + m->r1.start,__func__,__FILE__,__LINE__);

  }

  if (m->b)
    pfree (m->b + m->r1.start,__func__,__FILE__,__LINE__);
  pfree (m,__func__,__FILE__,__LINE__);
}
#line 230 "/home/lisergey/basilisk/src/grid/memindex/range.h"
void mem_assign (struct _Memindex * m, int i, int j, int k, int len, void * b)
{
  Memalloc mem[] = {{(void **)&m->b, sizeof(char ***)},
      {(void **)&m->r2, sizeof(Memrange)},
      {(void **)&m->r3, sizeof(Memrange *)},
      {NULL}};
  memrange_alloc (&m->r1, mem, i);
  Memalloc mem1[] = {{(void **)&m->b[i], sizeof(char **)},
       {(void **)&m->r3[i], sizeof(Memrange)},
       {NULL}};
  memrange_alloc (&m->r2[i], mem1, j);
  Memalloc mem2[] = {{(void **)&m->b[i][j], sizeof(char *)},
       {NULL}};
  memrange_alloc (&m->r3[i][j], mem2, k);
  ((m)->b[i][j][k]) = b;
}
#line 271 "/home/lisergey/basilisk/src/grid/memindex/range.h"
void mem_free (struct _Memindex * m, int i, int j, int k, int len)
{
  Memalloc mem[] = {{(void **)&m->b[i][j], sizeof(char *)},
      {NULL}};
  if (memrange_free (&m->r3[i][j], mem, k)) {
    Memalloc mem[] = {{(void **)&m->b[i], sizeof(char **)},
        {(void **)&m->r3[i], sizeof(Memrange)},
        {NULL}};
    if (memrange_free (&m->r2[i], mem, j)) {
      Memalloc mem[] = {{(void **)&m->b, sizeof(char ***)},
   {(void **)&m->r2, sizeof(Memrange)},
   {(void **)&m->r3, sizeof(Memrange *)},
   {NULL}};
      memrange_free (&m->r1, mem, i);
    }
  }
}
#line 318 "/home/lisergey/basilisk/src/grid/memindex/range.h"
#define foreach_mem(_m, _len, _i) {\
  Point point = {0};\
  for (point.i = max(Period.x*2, (_m)->r1.start);\
       point.i < min(_len - Period.x*2, (_m)->r1.end);\
       point.i += _i)\
    if ((_m)->b[point.i])\
      for (point.j = max(Period.y*2, (_m)->r2[point.i].start);\
    point.j < min(_len - Period.y*2, (_m)->r2[point.i].end);\
    point.j += _i)\
 if ((_m)->b[point.i][point.j])\
   for (point.k = max(Period.z*2, (_m)->r3[point.i][point.j].start);\
        point.k < min(_len - Period.z*2, (_m)->r3[point.i][point.j].end);\
        point.k += _i)\
     if ((_m)->b[point.i][point.j][point.k]) {\

#line 332

#define end_foreach_mem() }}
#line 7 "/home/lisergey/basilisk/src/grid/tree.h"
#line 24 "/home/lisergey/basilisk/src/grid/tree.h"
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


  , face_z = 1 << 2

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


  int k;

} IndexLevel;

typedef struct {
  IndexLevel * p;
  int n, nm;
} CacheLevel;

typedef struct {
  int i;

  int j;


  int k;

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






  return cube(_size(depth))*size;

}

static Layer * new_layer (int depth)
{
  Layer * l = ((Layer *) pmalloc ((1)*sizeof(Layer),__func__,__FILE__,__LINE__));
  l->len = _size (depth);
  if (depth == 0)
    l->pool = NULL;
  else {
    size_t size = sizeof(Cell) + datasize;


    l->pool = mempool_new (poolsize (depth, size), (1 << 3)*size);
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
  pfree (l,__func__,__FILE__,__LINE__);
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


  int k;

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
    c->p = (IndexLevel *) prealloc (c->p, (c->nm)*sizeof(IndexLevel),__func__,__FILE__,__LINE__);
  }
  c->p[c->n].i = p.i;

  c->p[c->n].j = p.j;


  c->p[c->n].k = p.k;

  c->n++;
}

static void cache_level_shrink (CacheLevel * c)
{
  if (c->nm > (c->n/128 + 1)*128) {
    c->nm = (c->n/128 + 1)*128;
    if (!(c->nm > c->n)) qassert ("/home/lisergey/basilisk/src/grid/tree.h", 200, "c->nm > c->n");
    c->p = (IndexLevel *) prealloc (c->p, sizeof (Index)*c->nm,__func__,__FILE__,__LINE__);
  }
}

static void cache_append (Cache * c, Point p, unsigned short flags)
{
  if (c->n >= c->nm) {
    c->nm += 128;
    c->p = (Index *) prealloc (c->p, (c->nm)*sizeof(Index),__func__,__FILE__,__LINE__);
  }
  c->p[c->n].i = p.i;

  c->p[c->n].j = p.j;


  c->p[c->n].k = p.k;

  c->p[c->n].level = p.level;
  c->p[c->n].flags = flags;
  c->n++;
}

void cache_shrink (Cache * c)
{
  cache_level_shrink ((CacheLevel *)c);
}
#line 261 "/home/lisergey/basilisk/src/grid/tree.h"
#define allocated(a,l,n) (((point.i+a) >= (((Tree *)grid)->L[point.level]->m)->r1.start && (point.i+a) < (((Tree *)grid)->L[point.level]->m->r1.end) && (((Tree *)grid)->L[point.level]->m)->b[point.i+a] && (point.j+l) >= (((Tree *)grid)->L[point.level]->m)->r2[point.i+a].start && (point.j+l) < (((Tree *)grid)->L[point.level]->m)->r2[point.i+a].end && (((Tree *)grid)->L[point.level]->m)->b[point.i+a][point.j+l] && (point.k+n) >= (((Tree *)grid)->L[point.level]->m)->r3[point.i+a][point.j+l].start && (point.k+n) < (((Tree *)grid)->L[point.level]->m)->r3[point.i+a][point.j+l].end && (((Tree *)grid)->L[point.level]->m)->b[point.i+a][point.j+l][point.k+n])\
                                          )\

#line 263

#define NEIGHBOR(a,l,n) (((((Tree *)grid)->L[point.level]->m)->b[point.i+a][point.j+l][point.k+n])\
                                       )\

#line 266

#define PARENT(a,l,n) (((((Tree *)grid)->L[point.level-1]->m)->b[(point.i+2)/2+a][(point.j+2)/2+l][(point.k+2)/2+n])\
\
\
                              )\

#line 271

#define allocated_child(a,l,n) (level < depth() &&\
         ((2*point.i-2 +a) >= (((Tree *)grid)->L[point.level+1]->m)->r1.start && (2*point.i-2 +a) < (((Tree *)grid)->L[point.level+1]->m->r1.end) && (((Tree *)grid)->L[point.level+1]->m)->b[2*point.i-2 +a] && (2*point.j-2 +l) >= (((Tree *)grid)->L[point.level+1]->m)->r2[2*point.i-2 +a].start && (2*point.j-2 +l) < (((Tree *)grid)->L[point.level+1]->m)->r2[2*point.i-2 +a].end && (((Tree *)grid)->L[point.level+1]->m)->b[2*point.i-2 +a][2*point.j-2 +l] && (2*point.k-2 +n) >= (((Tree *)grid)->L[point.level+1]->m)->r3[2*point.i-2 +a][2*point.j-2 +l].start && (2*point.k-2 +n) < (((Tree *)grid)->L[point.level+1]->m)->r3[2*point.i-2 +a][2*point.j-2 +l].end && (((Tree *)grid)->L[point.level+1]->m)->b[2*point.i-2 +a][2*point.j-2 +l][2*point.k-2 +n])\
\
\
                             )\

#line 277

#define CHILD(a,l,n) (((((Tree *)grid)->L[point.level+1]->m)->b[2*point.i-2 +a][2*point.j-2 +l][2*point.k-2 +n])\
\
\
                            )\

#line 282


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
    point.k + n,\
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
\
\
\
\
  struct { int x, y, z; } child = {\
    2*((point.i+2)%2)-1, 2*((point.j+2)%2)-1, 2*((point.k+2)%2)-1\
  };\
\
  NOT_UNUSED(child);\
  Point parent = point; NOT_UNUSED(parent);\
  parent.level--;\
  parent.i = (point.i + 2)/2;\
\
  parent.j = (point.j + 2)/2;\
\
\
  parent.k = (point.k + 2)/2;\
\
\
\
\
\
\

#line 341


#line 1 "grid/foreach_cell.h"
#line 1 "/home/lisergey/basilisk/src/grid/foreach_cell.h"
#line 66 "/home/lisergey/basilisk/src/grid/foreach_cell.h"
#define foreach_cell_root(root)\
  {\
    int ig = 0, jg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg);\
    Point point = {0};\
\
\
\
\
\
    int kg = 0; NOT_UNUSED(kg);\
    struct { int l, i, j, k, stage; } stack[20];\
\
    int _s = -1;\
    { _s++; stack[_s].l = 0; stack[_s].i = root.i; stack[_s].j = root.j; stack[_s].k = root.k; stack[_s].stage = 0; };\
    while (_s >= 0) {\
      int stage;\
      { point.level = stack[_s].l; point.i = stack[_s].i; point.j = stack[_s].j; point.k = stack[_s].k; stage = stack[_s].stage; _s--; };\
      if (!allocated (0,0,0))\
 continue;\
      switch (stage) {\
      case 0: {\
 POINT_VARIABLES;\
\

#line 89

#define end_foreach_cell_root()\
        if (point.level < grid->depth) {\
   { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 1; };\
          { _s++; stack[_s].l = point.level + 1; stack[_s].i = (2*point.i - 2); stack[_s].j = (2*point.j - 2); stack[_s].k = (2*point.k - 2); stack[_s].stage = 0; };\
        }\
        break;\
      }\
\
      case 1: { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 2; };\
       { _s++; stack[_s].l = point.level + 1; stack[_s].i = (2*point.i - 2); stack[_s].j = (2*point.j - 2); stack[_s].k = ((2*point.k - 2) + 1); stack[_s].stage = 0; }; break;\
      case 2: { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 3; };\
       { _s++; stack[_s].l = point.level + 1; stack[_s].i = (2*point.i - 2); stack[_s].j = ((2*point.j - 2) + 1); stack[_s].k = (2*point.k - 2); stack[_s].stage = 0; }; break;\
      case 3: { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 4; };\
       { _s++; stack[_s].l = point.level + 1; stack[_s].i = (2*point.i - 2); stack[_s].j = ((2*point.j - 2) + 1); stack[_s].k = ((2*point.k - 2) + 1); stack[_s].stage = 0; }; break;\
      case 4: { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 5; };\
       { _s++; stack[_s].l = point.level + 1; stack[_s].i = ((2*point.i - 2) + 1); stack[_s].j = (2*point.j - 2); stack[_s].k = (2*point.k - 2); stack[_s].stage = 0; }; break;\
      case 5: { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 6; };\
       { _s++; stack[_s].l = point.level + 1; stack[_s].i = ((2*point.i - 2) + 1); stack[_s].j = (2*point.j - 2); stack[_s].k = ((2*point.k - 2) + 1); stack[_s].stage = 0; }; break;\
      case 6: { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 7; };\
       { _s++; stack[_s].l = point.level + 1; stack[_s].i = ((2*point.i - 2) + 1); stack[_s].j = ((2*point.j - 2) + 1); stack[_s].k = (2*point.k - 2); stack[_s].stage = 0; }; break;\
      case 7: { _s++; stack[_s].l = point.level + 1; stack[_s].i = ((2*point.i - 2) + 1); stack[_s].j = ((2*point.j - 2) + 1); stack[_s].k = ((2*point.k - 2) + 1); stack[_s].stage = 0; }; break;\
\
      }\
    }\
  }\

#line 123


#define foreach_cell() {\
\
\
\
\
\
  Point root = {2,2,2,0};\
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
      for (root.k = 2*Period.z; root.k <= 2*(2 - Period.z); root.k++)\
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
\
\
    int kg = 0; NOT_UNUSED(kg);\
    struct { int l, i, j, k, stage; } stack[20];\
\
    int _s = -1;\
    { _s++; stack[_s].l = 0; stack[_s].i = root.i; stack[_s].j = root.j; stack[_s].k = root.k; stack[_s].stage = 0; };\
    while (_s >= 0) {\
      int stage;\
      { point.level = stack[_s].l; point.i = stack[_s].i; point.j = stack[_s].j; point.k = stack[_s].k; stage = stack[_s].stage; _s--; };\
      if (!allocated (0,0,0))\
 continue;\
      switch (stage) {\
      case 0: {\
        POINT_VARIABLES;\
 if (point.level == grid->depth) {\
   { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 8; };\
 }\
 else {\
   { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 1; };\
   if (condition)\
     { _s++; stack[_s].l = point.level + 1; stack[_s].i = (2*point.i - 2); stack[_s].j = (2*point.j - 2); stack[_s].k = (2*point.k - 2); stack[_s].stage = 0; };\
 }\
 break;\
      }\
\
      case 1:\
 { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 2; };\
 if (condition)\
   { _s++; stack[_s].l = point.level + 1; stack[_s].i = (2*point.i - 2); stack[_s].j = (2*point.j - 2); stack[_s].k = ((2*point.k - 2) + 1); stack[_s].stage = 0; };\
 break;\
      case 2:\
 { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 3; };\
 if (condition)\
   { _s++; stack[_s].l = point.level + 1; stack[_s].i = (2*point.i - 2); stack[_s].j = ((2*point.j - 2) + 1); stack[_s].k = (2*point.k - 2); stack[_s].stage = 0; };\
 break;\
      case 3:\
 { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 4; };\
 if (condition)\
   { _s++; stack[_s].l = point.level + 1; stack[_s].i = (2*point.i - 2); stack[_s].j = ((2*point.j - 2) + 1); stack[_s].k = ((2*point.k - 2) + 1); stack[_s].stage = 0; };\
 break;\
      case 4:\
 { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 5; };\
 if (condition)\
   { _s++; stack[_s].l = point.level + 1; stack[_s].i = ((2*point.i - 2) + 1); stack[_s].j = (2*point.j - 2); stack[_s].k = (2*point.k - 2); stack[_s].stage = 0; };\
 break;\
      case 5:\
 { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 6; };\
 if (condition)\
   { _s++; stack[_s].l = point.level + 1; stack[_s].i = ((2*point.i - 2) + 1); stack[_s].j = (2*point.j - 2); stack[_s].k = ((2*point.k - 2) + 1); stack[_s].stage = 0; };\
 break;\
      case 6:\
 { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 7; };\
 if (condition)\
   { _s++; stack[_s].l = point.level + 1; stack[_s].i = ((2*point.i - 2) + 1); stack[_s].j = ((2*point.j - 2) + 1); stack[_s].k = (2*point.k - 2); stack[_s].stage = 0; };\
 break;\
      case 7:\
 { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 8; };\
 if (condition)\
   { _s++; stack[_s].l = point.level + 1; stack[_s].i = ((2*point.i - 2) + 1); stack[_s].j = ((2*point.j - 2) + 1); stack[_s].k = ((2*point.k - 2) + 1); stack[_s].stage = 0; };\
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
\
\
    Point root = {2,2,2,0};\
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
      for (root.k = 0; root.k <= 2*2; root.k++)\
\
 foreach_cell_post_root (condition, root)\

#line 275

#define end_foreach_cell_post_all() end_foreach_cell_post_root() }

#define foreach_leaf() foreach_cell()\
  if (is_leaf (cell)) {\
    if (is_active(cell) && is_local(cell)) {\

#line 281

#define end_foreach_leaf() } continue; } end_foreach_cell()
#line 344 "/home/lisergey/basilisk/src/grid/tree.h"
#line 379 "/home/lisergey/basilisk/src/grid/tree.h"
#define foreach_child() {\
  int _i = 2*point.i - 2, _j = 2*point.j - 2, _k = 2*point.k - 2;\
  point.level++;\
  for (int _l = 0; _l < 2; _l++) {\
    point.i = _i + _l;\
    for (int _m = 0; _m < 2; _m++) {\
      point.j = _j + _m;\
      for (int _n = 0; _n < 2; _n++) {\
 point.k = _k + _n;\
 POINT_VARIABLES;\

#line 389

#define end_foreach_child()\
      }\
    }\
  }\
  point.i = (_i + 2)/2;point.j = (_j + 2)/2;point.k = (_k + 2)/2;\
  point.level--;\
}\

#line 397

#define foreach_child_break() _l = _m = _n = 2
#line 407 "/home/lisergey/basilisk/src/grid/tree.h"
#define is_refined_check() ((!is_leaf (cell) && cell.neighbors && cell.pid >= 0) &&\
    point.i > 0 && point.i < (1 << level) + 2*2 - 1\
\
    && point.j > 0 && point.j < (1 << level) + 2*2 - 1\
\
\
    && point.k > 0 && point.k < (1 << level) + 2*2 - 1\
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
  point.k = 2;\
\
  int _k; unsigned short _flags; NOT_UNUSED(_flags);\
  OMP(omp for schedule(static))\
  for (_k = 0; _k < _cache.n; _k++) {\
    point.i = _cache.p[_k].i;\
\
    point.j = _cache.p[_k].j;\
\
\
    point.k = _cache.p[_k].k;\
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
  point.k = 2;\
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
    point.k = _cache.p[_k].k;\
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
 for (int _d = 0; _d < 3; _d++) {\
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
       point.k -= kg; z -= kg*Delta/2.;\
\

#line 499

#define end_foreach_boundary()\
       point.i += ig; x += ig*Delta/2.;\
\
       point.j += jg; y += jg*Delta/2.;\
\
\
       point.k += kg; z += kg*Delta/2.;\
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
#line 1 "/home/lisergey/basilisk/src/grid/neighbors.h"
#line 35 "/home/lisergey/basilisk/src/grid/neighbors.h"
#define foreach_neighbor(_s) {\
  int _nn = _s + 0 ? _s + 0 : 2;\
  int _i = point.i, _j = point.j, _k = point.k;\
  for (int _l = - _nn; _l <= _nn; _l++) {\
    point.i = _i + _l;\
    for (int _m = - _nn; _m <= _nn; _m++) {\
      point.j = _j + _m;\
      for (int _n = - _nn; _n <= _nn; _n++) {\
 point.k = _k + _n;\
 POINT_VARIABLES;\

#line 45

#define end_foreach_neighbor()\
      }\
    }\
  }\
  point.i = _i; point.j = _j; point.k = _k;\
}\

#line 52

#define foreach_neighbor_break() _l = _m = _n = _nn + 1
#line 524 "/home/lisergey/basilisk/src/grid/tree.h"

static inline bool has_local_children (Point point)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;
  {foreach_child()
    if (is_local(cell))
      return true;end_foreach_child()}
  return false;
}

static inline void cache_append_face (Point point, unsigned short flags)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;
  Tree * q = ((Tree *)grid);
  cache_append (&q->faces, point, flags);
#line 548 "/home/lisergey/basilisk/src/grid/tree.h"
  
    if (flags & face_x)
      for (int i = 0; i <= 1; i++)
 for (int j = 0; j <= 1; j++)
   if (!is_vertex(neighbor(0,i,j))) {
     cache_append (&q->vertices, neighborp(0,i,j), 0);
     neighbor(0,i,j).flags |= vertex;
   }
    
#line 549
if (flags & face_y)
      for (int i = 0; i <= 1; i++)
 for (int j = 0; j <= 1; j++)
   if (!is_vertex(neighbor(j,0,i))) {
     cache_append (&q->vertices, neighborp(j,0,i), 0);
     neighbor(j,0,i).flags |= vertex;
   }
    
#line 549
if (flags & face_z)
      for (int i = 0; i <= 1; i++)
 for (int j = 0; j <= 1; j++)
   if (!is_vertex(neighbor(i,j,0))) {
     cache_append (&q->vertices, neighborp(i,j,0), 0);
     neighbor(i,j,0).flags |= vertex;
   }

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
#line 601 "/home/lisergey/basilisk/src/grid/tree.h"
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
   
#line 619
if ((neighbor(0,0,-1).pid < 0) || (!is_leaf(neighbor(0,0,-1)) && !neighbor(0,0,-1).neighbors && neighbor(0,0,-1).pid >= 0) ||
       is_leaf(neighbor(0,0,-1)))
     flags |= face_z;
 if (flags)
   cache_append (&q->faces, point, flags);
 
   if ((neighbor(1,0,0).pid < 0) || (!is_leaf(neighbor(1,0,0)) && !neighbor(1,0,0).neighbors && neighbor(1,0,0).pid >= 0) ||
       (!is_local(neighbor(1,0,0)) && is_leaf(neighbor(1,0,0))))
     cache_append (&q->faces, neighborp(1,0,0), face_x);
   
#line 625
if ((neighbor(0,1,0).pid < 0) || (!is_leaf(neighbor(0,1,0)) && !neighbor(0,1,0).neighbors && neighbor(0,1,0).pid >= 0) ||
       (!is_local(neighbor(0,1,0)) && is_leaf(neighbor(0,1,0))))
     cache_append (&q->faces, neighborp(0,1,0), face_y);
   
#line 625
if ((neighbor(0,0,1).pid < 0) || (!is_leaf(neighbor(0,0,1)) && !neighbor(0,0,1).neighbors && neighbor(0,0,1).pid >= 0) ||
       (!is_local(neighbor(0,0,1)) && is_leaf(neighbor(0,0,1))))
     cache_append (&q->faces, neighborp(0,0,1), face_z);

 for (int i = 0; i <= 1; i++)

   for (int j = 0; j <= 1; j++)


     for (int k = 0; k <= 1; k++)

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
   
#line 648
if (allocated(0,0,-1) &&
       is_local(neighbor(0,0,-1)) && (!is_leaf(neighbor(0,0,-1)) && !neighbor(0,0,-1).neighbors && neighbor(0,0,-1).pid >= 0))
     flags |= face_z;
 if (flags)
   cache_append_face (point, flags);
 
   if (allocated(1,0,0) && is_local(neighbor(1,0,0)) &&
       (!is_leaf(neighbor(1,0,0)) && !neighbor(1,0,0).neighbors && neighbor(1,0,0).pid >= 0))
     cache_append_face (neighborp(1,0,0), face_x);
   
#line 654
if (allocated(0,1,0) && is_local(neighbor(0,1,0)) &&
       (!is_leaf(neighbor(0,1,0)) && !neighbor(0,1,0).neighbors && neighbor(0,1,0).pid >= 0))
     cache_append_face (neighborp(0,1,0), face_y);
   
#line 654
if (allocated(0,0,1) && is_local(neighbor(0,0,1)) &&
       (!is_leaf(neighbor(0,0,1)) && !neighbor(0,0,1).neighbors && neighbor(0,0,1).pid >= 0))
     cache_append_face (neighborp(0,0,1), face_z);
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

#if !1
  grid->tn = grid->n;
  grid->maxdepth = grid->depth;
#endif
}

#define foreach() { if (((Tree *)grid)->dirty) update_cache_f(); }; foreach_cache(((Tree *)grid)->leaves)
#define end_foreach() end_foreach_cache()

#define foreach_face_generic()\
  { if (((Tree *)grid)->dirty) update_cache_f(); };\
  foreach_cache(((Tree *)grid)->faces) 
#line 736

#define end_foreach_face_generic() end_foreach_cache()

#define is_face_x() { int ig = -1; VARIABLES; if (_flags & face_x) {
#define end_is_face_x() }}


#define is_face_y() { int jg = -1; VARIABLES; if (_flags & face_y) {
#define end_is_face_y() }}


#define is_face_z() { int kg = -1; VARIABLES; if (_flags & face_z) {
#define end_is_face_z() }}


#define foreach_vertex()\
  { if (((Tree *)grid)->dirty) update_cache_f(); };\
  foreach_cache(((Tree *)grid)->vertices) {\
    x -= Delta/2.;\
\
    y -= Delta/2.;\
\
\
    z -= Delta/2.;\
\

#line 762

#define end_foreach_vertex() } end_foreach_cache()
#line 734 "/home/lisergey/basilisk/src/grid/tree.h"
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
    pfree (name[i].p,__func__,__FILE__,__LINE__);
  pfree (name,__func__,__FILE__,__LINE__);
  return ((CacheLevel *) pcalloc (depth() + 1, sizeof(CacheLevel),__func__,__FILE__,__LINE__));
}

static void update_depth (int inc)
{
  Tree * q = ((Tree *)grid);
  grid->depth += inc;
  q->L = &(q->L[-1]);
  q->L = (Layer * *) prealloc (q->L, (grid->depth + 2)*sizeof(Layer *),__func__,__FILE__,__LINE__);
  q->L = &(q->L[1]);
  if (inc > 0)
    q->L[grid->depth] = new_layer (grid->depth);
  q->active = cache_level_resize (q->active, inc);
  q->prolongation = cache_level_resize (q->prolongation, inc);
  q->boundary = cache_level_resize (q->boundary, inc);
  q->restriction = cache_level_resize (q->restriction, inc);
}
#line 861 "/home/lisergey/basilisk/src/grid/tree.h"
typedef void (* PeriodicFunction) (struct _Memindex *, int, int, int, int, void *);

static void periodic_function (struct _Memindex * m, int i, int j, int k, int len,
          void * b, PeriodicFunction f)
{
  f(m, i, j, k, len, b);
  if (Period.x) {
    int nl = len - 2*2;
    for (int l = - 1; l <= 1; l += 2)
      for (int n = i + l*nl; n >= 0 && n < len; n += l*nl)
 f(m, n, j, k, len, b);
    if (Period.y) {
      for (int l = - 1; l <= 1; l += 2)
 for (int n = j + l*nl; n >= 0 && n < len; n += l*nl) {
   f(m, i, n, k, len, b);
   for (int o = - 1; o <= 1; o += 2)
     for (int p = i + o*nl; p >= 0 && p < len; p += o*nl)
       f(m, p, n, k, len, b);
 }
      if (Period.z)
 for (int l = - 1; l <= 1; l += 2)
   for (int n = k + l*nl; n >= 0 && n < len; n += l*nl) {
     f(m, i, j, n, len, b);
     for (int q = - 1; q <= 1; q += 2)
       for (int r = j + q*nl; r >= 0 && r < len; r += q*nl)
  f(m, i, r, n, len, b);
     for (int o = - 1; o <= 1; o += 2)
       for (int p = i + o*nl; p >= 0 && p < len; p += o*nl) {
  f(m, p, j, n, len, b);
  for (int q = - 1; q <= 1; q += 2)
    for (int r = j + q*nl; r >= 0 && r < len; r += q*nl)
      f(m, p, r, n, len, b);
       }
   }
    }
    else if (Period.z)
      for (int l = - 1; l <= 1; l += 2)
 for (int n = k + l*nl; n >= 0 && n < len; n += l*nl) {
   f(m, i, j, n, len, b);
   for (int o = - 1; o <= 1; o += 2)
     for (int p = i + o*nl; p >= 0 && p < len; p += o*nl)
       f(m, p, j, n, len, b);
 }
  }
  else if (Period.y) {
    int nl = len - 2*2;
    for (int l = - 1; l <= 1; l += 2)
      for (int n = j + l*nl; n >= 0 && n < len; n += l*nl)
 f(m, i, n, k, len, b);
    if (Period.z)
      for (int l = - 1; l <= 1; l += 2)
 for (int n = k + l*nl; n >= 0 && n < len; n += l*nl) {
   f(m, i, j, n, len, b);
   for (int o = - 1; o <= 1; o += 2)
     for (int p = j + o*nl; p >= 0 && p < len; p += o*nl)
       f(m, i, p, n, len, b);
 }
  }
  else if (Period.z) {
    int nl = len - 2*2;
    for (int l = - 1; l <= 1; l += 2)
      for (int n = k + l*nl; n >= 0 && n < len; n += l*nl)
 f(m, i, j, n, len, b);
  }
}

static void assign_periodic (struct _Memindex * m, int i, int j, int k, int len, void * b)
{
  periodic_function (m, i, j, k, len, b, mem_assign);
}

static void free_periodic (struct _Memindex * m, int i, int j, int k, int len)
{
  periodic_function (m, i, j, k, len, NULL, (PeriodicFunction) mem_free);
}


static void alloc_children (Point point)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;
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
#line 962 "/home/lisergey/basilisk/src/grid/tree.h"
    int j = 2*point.j - 2;
    for (int l = 0; l < 2; l++, j++) {
      int m = 2*point.k - 2;
      for (int n = 0; n < 2; n++, m++) {
 assign_periodic (L->m, i, j, m, L->len, b);
 b += len;
      }
    }

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
#line 1017 "/home/lisergey/basilisk/src/grid/tree.h"
static void free_children (Point point)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;

  Layer * L = ((Tree *)grid)->L[point.level + 1];
  int i = 2*point.i - 2;
  if (!(((L->m)->b[i][2*point.j - 2][2*point.k - 2]))) qassert ("/home/lisergey/basilisk/src/grid/tree.h", 1022, "mem_data (L->m,i,2*point.j - GHOSTS,2*point.k - GHOSTS)");
  mempool_free (L->pool, ((L->m)->b[i][2*point.j - 2][2*point.k - 2])
                                               );
  for (int k = 0; k < 2; k++, i++) {
    int j = 2*point.j - 2;
    for (int l = 0; l < 2; l++, j++) {
      int m = 2*point.k - 2;
      for (int n = 0; n < 2; n++, m++)
 free_periodic (L->m, i, j, m, L->len);
    }
  }
  if (--L->nc == 0) {
    destroy_layer (L);
    if (!(point.level + 1 == grid->depth)) qassert ("/home/lisergey/basilisk/src/grid/tree.h", 1035, "point.level + 1 == grid->depth");
    update_depth (-1);
  }
}


void increment_neighbors (Point point)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;
  ((Tree *)grid)->dirty = true;
  if (cell.neighbors++ == 0)
    alloc_children (point);
  {foreach_neighbor (2/2)
    if (cell.neighbors++ == 0)
      alloc_children (point);end_foreach_neighbor()}
  cell.neighbors--;
}

void decrement_neighbors (Point point)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;
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
#line 1088 "/home/lisergey/basilisk/src/grid/tree.h"
    char * p = (char *) prealloc (((L->m)->b[point.i][point.j][point.k]),
     newlen*sizeof(char),__func__,__FILE__,__LINE__);
    assign_periodic (L->m, point.i, point.j, point.k, L->len, p);

  }end_foreach_mem();}

  for (int l = 1; l <= depth(); l++) {
    Layer * L = q->L[l];
    Mempool * oldpool = L->pool;
    L->pool = mempool_new (poolsize (l, newlen), (1 << 3)*newlen);
    {foreach_mem (L->m, L->len, 2) {
      char * new = (char *) mempool_alloc (L->pool);
#line 1114 "/home/lisergey/basilisk/src/grid/tree.h"
      for (int l = 0; l < 2; l++)
 for (int m = 0; m < 2; m++)
   for (int n = 0; n < 2; n++) {
     memcpy (new, ((L->m)->b[point.i + l][point.j + m][point.k + n]),
      oldlen);
     assign_periodic (L->m, point.i + l, point.j + m, point.k + n,
     L->len, new);
     new += newlen;
   }

    }end_foreach_mem();}
    mempool_destroy (oldpool);
  }
}



#define VN v.x
#define VT v.y
#define VR v.z




#if 1
# define disable_fpe_for_mpi() disable_fpe (FE_DIVBYZERO|FE_INVALID)
# define enable_fpe_for_mpi() enable_fpe (FE_DIVBYZERO|FE_INVALID)
#else
# define disable_fpe_for_mpi()
# define enable_fpe_for_mpi()
#endif

static inline void no_restriction (Point point, scalar s);

static bool normal_neighbor (Point point, scalar * scalars, vector * vectors)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;
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


       scalar vr = VR;
       val(v.z,0,0,0) = _attribute[vr.i].boundary[id](neighbor, point, v.z, NULL);

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
       val(v.z,0,0,0) = _attribute[vt.i].boundary[id](neighbor, point, v.z, NULL);


       scalar vr = VR;
       val(v.x,0,0,0) = _attribute[vr.i].boundary[id](neighbor, point, v.x, NULL);

     }}}
   return true;
 }
      
#line 1152
for (int i = -k; i <= k; i += 2*k)
 if ((allocated(0,0,i) && !(neighbor(0,0,i).pid < 0))) {
   Point neighbor = neighborp(0,0,i);
   int id = (- cell.pid - 1);
   {scalar*_i=(scalar*)( scalars);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    
       val(s,0,0,0) = _attribute[s.i].boundary[id](neighbor, point, s, NULL);}}
   {vector*_i=(vector*)( vectors);if(_i)for(vector v=*_i;(&v)->x.i>=0;v=*++_i){
     {
       scalar vn = VN;
       val(v.z,0,0,0) = _attribute[vn.i].boundary[id](neighbor, point, v.z, NULL);

       scalar vt = VT;
       val(v.x,0,0,0) = _attribute[vt.i].boundary[id](neighbor, point, v.x, NULL);


       scalar vr = VR;
       val(v.y,0,0,0) = _attribute[vr.i].boundary[id](neighbor, point, v.y, NULL);

     }}}
   return true;
 }}
  return false;
}

static bool diagonal_neighbor_2D (Point point,
      scalar * scalars, vector * vectors)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;

  for (int k = 1; k <= 2; k++)

    {

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

  scalar vr = VR;
  val(v.z,0,0,0) = (_attribute[vr.i].boundary[id1](n,n1,v.z,NULL) +
    _attribute[vr.i].boundary[id2](n,n2,v.z,NULL) -
    val(v.z,i,j,0));

       }}}
     return true;
   }

      
#line 1185
for (int i = -k; i <= k; i += 2*k)
 for (int j = -k; j <= k; j += 2*k)
   if (allocated(0,i,j) && (allocated(0,i,j) && !(neighbor(0,i,j).pid < 0)) &&
       allocated(0,i,0) && (neighbor(0,i,0).pid < 0) &&
       allocated(0,0,j) && (neighbor(0,0,j).pid < 0)) {
     Point n = neighborp(0,i,j),
       n1 = neighborp(0,i,0), n2 = neighborp(0,0,j);
     int id1 = (- neighbor(0,i,0).pid - 1), id2 = (- neighbor(0,0,j).pid - 1);
     {scalar*_i=(scalar*)( scalars);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      
  val(s,0,0,0) = (_attribute[s.i].boundary[id1](n,n1,s,NULL) +
         _attribute[s.i].boundary[id2](n,n2,s,NULL) -
         val(s,0,i,j));}}
     {vector*_i=(vector*)( vectors);if(_i)for(vector v=*_i;(&v)->x.i>=0;v=*++_i){
       {
  scalar vt = VT, vn = VN;
  val(v.y,0,0,0) = (_attribute[vt.i].boundary[id1](n,n1,v.y,NULL) +
    _attribute[vn.i].boundary[id2](n,n2,v.y,NULL) -
    val(v.y,0,i,j));
  val(v.z,0,0,0) = (_attribute[vn.i].boundary[id1](n,n1,v.z,NULL) +
    _attribute[vt.i].boundary[id2](n,n2,v.z,NULL) -
    val(v.z,0,i,j));

  scalar vr = VR;
  val(v.x,0,0,0) = (_attribute[vr.i].boundary[id1](n,n1,v.x,NULL) +
    _attribute[vr.i].boundary[id2](n,n2,v.x,NULL) -
    val(v.x,0,i,j));

       }}}
     return true;
   }

      
#line 1185
for (int i = -k; i <= k; i += 2*k)
 for (int j = -k; j <= k; j += 2*k)
   if (allocated(j,0,i) && (allocated(j,0,i) && !(neighbor(j,0,i).pid < 0)) &&
       allocated(0,0,i) && (neighbor(0,0,i).pid < 0) &&
       allocated(j,0,0) && (neighbor(j,0,0).pid < 0)) {
     Point n = neighborp(j,0,i),
       n1 = neighborp(0,0,i), n2 = neighborp(j,0,0);
     int id1 = (- neighbor(0,0,i).pid - 1), id2 = (- neighbor(j,0,0).pid - 1);
     {scalar*_i=(scalar*)( scalars);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      
  val(s,0,0,0) = (_attribute[s.i].boundary[id1](n,n1,s,NULL) +
         _attribute[s.i].boundary[id2](n,n2,s,NULL) -
         val(s,j,0,i));}}
     {vector*_i=(vector*)( vectors);if(_i)for(vector v=*_i;(&v)->x.i>=0;v=*++_i){
       {
  scalar vt = VT, vn = VN;
  val(v.z,0,0,0) = (_attribute[vt.i].boundary[id1](n,n1,v.z,NULL) +
    _attribute[vn.i].boundary[id2](n,n2,v.z,NULL) -
    val(v.z,j,0,i));
  val(v.x,0,0,0) = (_attribute[vn.i].boundary[id1](n,n1,v.x,NULL) +
    _attribute[vt.i].boundary[id2](n,n2,v.x,NULL) -
    val(v.x,j,0,i));

  scalar vr = VR;
  val(v.y,0,0,0) = (_attribute[vr.i].boundary[id1](n,n1,v.y,NULL) +
    _attribute[vr.i].boundary[id2](n,n2,v.y,NULL) -
    val(v.y,j,0,i));

       }}}
     return true;
   }}

  return false;
}

static bool diagonal_neighbor_3D (Point point,
      scalar * scalars, vector * vectors)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;

  for (int n = 1; n <= 2; n++)
    for (int i = -n; i <= n; i += 2*n)
      for (int j = -n; j <= n; j += 2*n)
 for (int k = -n; k <= n; k += 2*n)
   if ((allocated(i,j,k) && !(neighbor(i,j,k).pid < 0)) &&
       (neighbor(i,j,0).pid < 0) &&
       (neighbor(i,0,k).pid < 0) &&
       (neighbor(0,j,k).pid < 0)) {
     Point
       n0 = neighborp(i,j,k),
       n1 = neighborp(i,j,0),
       n2 = neighborp(i,0,k),
       n3 = neighborp(0,j,k);
     int
       id1 = (- neighbor(i,j,0).pid - 1),
       id2 = (- neighbor(i,0,k).pid - 1),
       id3 = (- neighbor(0,j,k).pid - 1);
     {scalar*_i=(scalar*)( scalars);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      
  val(s,0,0,0) = (_attribute[s.i].boundary[id1](n0,n1,s,NULL) +
         _attribute[s.i].boundary[id2](n0,n2,s,NULL) +
         _attribute[s.i].boundary[id3](n0,n3,s,NULL) -
         2.*val(s,i,j,k));}}
     {vector*_i=(vector*)( vectors);if(_i)for(vector v=*_i;(&v)->x.i>=0;v=*++_i){
       {
  scalar vt = VT, vn = VN, vr = VR;
  val(v.x,0,0,0) = (_attribute[vt.i].boundary[id1](n0,n1,v.x,NULL) +
    _attribute[vt.i].boundary[id2](n0,n2,v.x,NULL) +
    _attribute[vn.i].boundary[id3](n0,n3,v.x,NULL) -
    2.*val(v.x,i,j,k));
  val(v.y,0,0,0) = (_attribute[vt.i].boundary[id1](n0,n1,v.y,NULL) +
    _attribute[vn.i].boundary[id2](n0,n2,v.y,NULL) +
    _attribute[vt.i].boundary[id3](n0,n3,v.y,NULL) -
    2.*val(v.y,i,j,k));
  val(v.z,0,0,0) = (_attribute[vn.i].boundary[id1](n0,n1,v.z,NULL) +
    _attribute[vr.i].boundary[id2](n0,n2,v.z,NULL) +
    _attribute[vr.i].boundary[id3](n0,n3,v.z,NULL) -
    2.*val(v.z,i,j,k));
       }}}
     return true;
   }

  return false;
}



static Point tangential_neighbor_x (Point point, bool * zn)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;
  for (int k = 1; k <= 2; k++)
    for (int j = -k; j <= k; j += 2*k) {
      if ((allocated(0,j,0) && !(neighbor(0,j,0).pid < 0)) || (allocated(-1,j,0) && !(neighbor(-1,j,0).pid < 0))) {
 *zn = false;
 return neighborp(0,j,0);
      }


      if ((allocated(0,0,j) && !(neighbor(0,0,j).pid < 0)) || (allocated(-1,0,j) && !(neighbor(-1,0,j).pid < 0))) {
 *zn = true;
 return neighborp(0,0,j);
      }

    }
  return (Point){.level = -1};
}

#line 1271
static Point tangential_neighbor_y (Point point, bool * zn)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;
  for (int k = 1; k <= 2; k++)
    for (int j = -k; j <= k; j += 2*k) {
      if ((allocated(0,0,j) && !(neighbor(0,0,j).pid < 0)) || (allocated(0,-1,j) && !(neighbor(0,-1,j).pid < 0))) {
 *zn = false;
 return neighborp(0,0,j);
      }


      if ((allocated(j,0,0) && !(neighbor(j,0,0).pid < 0)) || (allocated(j,-1,0) && !(neighbor(j,-1,0).pid < 0))) {
 *zn = true;
 return neighborp(j,0,0);
      }

    }
  return (Point){.level = -1};
}

#line 1271
static Point tangential_neighbor_z (Point point, bool * zn)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;
  for (int k = 1; k <= 2; k++)
    for (int j = -k; j <= k; j += 2*k) {
      if ((allocated(j,0,0) && !(neighbor(j,0,0).pid < 0)) || (allocated(j,0,-1) && !(neighbor(j,0,-1).pid < 0))) {
 *zn = false;
 return neighborp(j,0,0);
      }


      if ((allocated(0,j,0) && !(neighbor(0,j,0).pid < 0)) || (allocated(0,j,-1) && !(neighbor(0,j,-1).pid < 0))) {
 *zn = true;
 return neighborp(0,j,0);
      }

    }
  return (Point){.level = -1};
}


static inline bool is_boundary_point (Point point) {int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;
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
val(v.y,0,0,0) = undefined;
     
#line 1323
val(v.z,0,0,0) = undefined;}}}
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



  scalar vt = zn ? VT : VR;

 
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



  scalar vt = zn ? VT : VR;

 
    val(v.y,0,0,0) = _attribute[vt.i].boundary[id](neighbor, point, v.y, NULL);
       }}}
     }
     else

       {vector*_i=(vector*)( faces);if(_i)for(vector v=*_i;(&v)->x.i>=0;v=*++_i){
 
    val(v.y,0,0,0) = 0.;}}
   }

 }
 
#line 1328
for (int i = -1; i <= 1; i += 2) {

   if ((allocated(0,0,i) && !(neighbor(0,0,i).pid < 0))) {
     Point neighbor = neighborp(0,0,i);
     {vector*_i=(vector*)( faces);if(_i)for(vector v=*_i;(&v)->x.i>=0;v=*++_i){ {
       scalar vn = VN;
       if (_attribute[vn.i].boundary[id])
 
    val(v.z,0,0,(i + 1)/2) = _attribute[vn.i].boundary[id](neighbor, point, v.z, NULL);
     }}}
   }

   else if (i == -1) {

     bool zn;
     Point neighbor = tangential_neighbor_z (point, &zn);
     if (neighbor.level >= 0) {
       int id = is_boundary_point (neighbor) ?
  (- neighbor(0,0,-1).pid - 1) : (- cell.pid - 1);
       {vector*_i=(vector*)( faces);if(_i)for(vector v=*_i;(&v)->x.i>=0;v=*++_i){ {



  scalar vt = zn ? VT : VR;

 
    val(v.z,0,0,0) = _attribute[vt.i].boundary[id](neighbor, point, v.z, NULL);
       }}}
     }
     else

       {vector*_i=(vector*)( faces);if(_i)for(vector v=*_i;(&v)->x.i>=0;v=*++_i){
 
    val(v.z,0,0,0) = 0.;}}
   }

 }
    }
  }end_foreach_boundary_level();}

  pfree (scalars,__func__,__FILE__,__LINE__);
  pfree (vectors,__func__,__FILE__,__LINE__);
  pfree (faces,__func__,__FILE__,__LINE__);
  enable_fpe_for_mpi();
}



#undef VN
#undef VT
#define VN _attribute[s.i].v.x
#define VT _attribute[s.i].v.y

static double masked_average (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;
  double sum = 0., n = 0.;
  {foreach_child()
    if (!(cell.pid < 0) && val(s,0,0,0) != 1e30)
      sum += val(s,0,0,0), n++;end_foreach_child()}
  return n ? sum/n : 1e30;
}


static double masked_average_x (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;
  double sum = 0., n = 0.;
  {foreach_child()
    if (child.x < 0 && (!(cell.pid < 0) || !(neighbor(1,0,0).pid < 0)) &&
 val(s,1,0,0) != 1e30)
      sum += val(s,1,0,0), n++;end_foreach_child()}
  return n ? sum/n : 1e30;
}

#line 1391
static double masked_average_y (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;
  double sum = 0., n = 0.;
  {foreach_child()
    if (child.y < 0 && (!(cell.pid < 0) || !(neighbor(0,1,0).pid < 0)) &&
 val(s,0,1,0) != 1e30)
      sum += val(s,0,1,0), n++;end_foreach_child()}
  return n ? sum/n : 1e30;
}

#line 1391
static double masked_average_z (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;
  double sum = 0., n = 0.;
  {foreach_child()
    if (child.z < 0 && (!(cell.pid < 0) || !(neighbor(0,0,1).pid < 0)) &&
 val(s,0,0,1) != 1e30)
      sum += val(s,0,0,1), n++;end_foreach_child()}
  return n ? sum/n : 1e30;
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
      } 
#line 1418
{
 double average = masked_average_z (parent, v.z);
 if ((neighbor(0,0,-1).pid < 0))
   val(v.z,0,0,0) = average;
 if ((neighbor(0,0,1).pid < 0))
   val(v.z,0,0,1) = average;
      }}}}
  }end_foreach_halo();}

  pfree (scalars,__func__,__FILE__,__LINE__);
  pfree (faces,__func__,__FILE__,__LINE__);
}
#line 1454 "/home/lisergey/basilisk/src/grid/tree.h"
static void free_cache (CacheLevel * c)
{
  for (int l = 0; l <= depth(); l++)
    pfree (c[l].p,__func__,__FILE__,__LINE__);
  pfree (c,__func__,__FILE__,__LINE__);
}

void free_grid (void)
{
  if (!grid)
    return;
  free_boundaries();
  Tree * q = ((Tree *)grid);
  pfree (q->leaves.p,__func__,__FILE__,__LINE__);
  pfree (q->faces.p,__func__,__FILE__,__LINE__);
  pfree (q->vertices.p,__func__,__FILE__,__LINE__);
  pfree (q->refined.p,__func__,__FILE__,__LINE__);


  Layer * L = q->L[0];
  {foreach_mem (L->m, L->len, 1) {





    pfree (((L->m)->b[point.i][point.j][point.k]),__func__,__FILE__,__LINE__);

  }end_foreach_mem();}
  for (int l = 0; l <= depth(); l++)
    destroy_layer (q->L[l]);
  q->L = &(q->L[-1]);
  pfree (q->L,__func__,__FILE__,__LINE__);
  free_cache (q->active);
  free_cache (q->prolongation);
  free_cache (q->boundary);
  free_cache (q->restriction);
  pfree (q,__func__,__FILE__,__LINE__);
  grid = NULL;
}

static void refine_level (int depth);

     
void init_grid (int n)
{tracing("init_grid","/home/lisergey/basilisk/src/grid/tree.h",1498);

  if (!(sizeof(Cell) % 8 == 0)) qassert ("/home/lisergey/basilisk/src/grid/tree.h", 1501, "sizeof(Cell) % 8 == 0");

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
  Tree * q = ((Tree *) pcalloc (1, sizeof(Tree),__func__,__FILE__,__LINE__));
  grid = (Grid *) q;
  grid->depth = 0;


  q->L = ((Layer * *) pmalloc ((2)*sizeof(Layer *),__func__,__FILE__,__LINE__));

  q->L[0] = NULL; q->L = &(q->L[1]);

  Layer * L = new_layer (0);
  q->L[0] = L;
#line 1554 "/home/lisergey/basilisk/src/grid/tree.h"
  for (int i = Period.x*2; i < L->len - Period.x*2; i++)
    for (int j = Period.y*2; j < L->len - Period.y*2; j++)
      for (int k = Period.z*2; k < L->len - Period.z*2; k++)
 assign_periodic (L->m, i, j, k, L->len,
    (char *) pcalloc (1, sizeof(Cell) + datasize,__func__,__FILE__,__LINE__));
  CELL(((L->m)->b[2][2][2])).flags |= leaf;
  if (pid() == 0)
    CELL(((L->m)->b[2][2][2])).flags |= active;
  for (int k = - 2*(1 - Period.x); k <= 2*(1 - Period.x); k++)
    for (int l = -2*(1 - Period.y); l <= 2*(1 - Period.y); l++)
      for (int n = -2*(1 - Period.z); n <= 2*(1 - Period.z); n++)
 CELL(((L->m)->b[2 +k][2 +l][2 +n])).pid =
   (k > 0 ? -1 - right :
    k < 0 ? -1 - left :
    l > 0 ? -1 - top :
    l < 0 ? -1 - bottom :
    n > 0 ? -1 - front :
    n < 0 ? -1 - back :
    0);
  CELL(((L->m)->b[2][2][2])).pid = 0;

  q->active = ((CacheLevel *) pcalloc (1, sizeof(CacheLevel),__func__,__FILE__,__LINE__));
  q->prolongation = ((CacheLevel *) pcalloc (1, sizeof(CacheLevel),__func__,__FILE__,__LINE__));
  q->boundary = ((CacheLevel *) pcalloc (1, sizeof(CacheLevel),__func__,__FILE__,__LINE__));
  q->restriction = ((CacheLevel *) pcalloc (1, sizeof(CacheLevel),__func__,__FILE__,__LINE__));
  q->dirty = true;
  N = 1 << depth;
#if 1
  void mpi_boundary_new();
  mpi_boundary_new();
#endif

  Boundary * b = ((Boundary *) pcalloc (1, sizeof(Boundary),__func__,__FILE__,__LINE__));
  b->level = box_boundary_level;
  b->restriction = masked_boundary_restriction;
  add_boundary (b);
  refine_level (depth);
  reset (all, 0.);
  { if (((Tree *)grid)->dirty) update_cache_f(); };
end_tracing("init_grid","/home/lisergey/basilisk/src/grid/tree.h",1593);}
#line 1629 "/home/lisergey/basilisk/src/grid/tree.h"
Point locate (double xp, double yp, double zp)
{
  for (int l = depth(); l >= 0; l--) {
    Point point = {0};int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;
    point.level = l;
    int n = 1 << point.level;
    point.i = (xp - X0)/L0*n + 2;

    point.j = (yp - Y0)/L0*n + 2;


    point.k = (zp - Z0)/L0*n + 2;

    if (point.i >= 0 && point.i < n + 2*2

 && point.j >= 0 && point.j < n + 2*2


 && point.k >= 0 && point.k < n + 2*2

 ) {
      if (allocated(0,0,0) && is_local(cell) && is_leaf(cell))
 return point;
    }
    else
      break;
  }
  Point point = {0};int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;
  point.level = -1;
  return point;
}



bool tree_is_full()
{
  { if (((Tree *)grid)->dirty) update_cache_f(); };
  return (grid->tn == 1L << grid->maxdepth*3);
}

#line 1 "grid/tree-common.h"
#line 1 "/home/lisergey/basilisk/src/grid/tree-common.h"



#line 1 "grid/multigrid-common.h"
#line 1 "/home/lisergey/basilisk/src/grid/multigrid-common.h"


#line 1 "grid/cartesian-common.h"
#line 1 "/home/lisergey/basilisk/src/grid/cartesian-common.h"
#line 1 "grid/events.h"
#line 1 "/home/lisergey/basilisk/src/grid/events.h"




static int END_EVENT = 1234567890;
static double TEND_EVENT = 1234567890;
static double TEPS = 1e-9;

static void event_error (Event * ev, const char * s)
{
  fprintf (ferr, "%s:%d: error: %s\n", ev->file, ev->line, s);
  exit (1);
}

static void init_event (Event * ev)
{
  if (ev->arrayi || ev->arrayt) {
    ev->i = -1; ev->t = - TEND_EVENT;
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
 int i = -123456; double t = - TEND_EVENT;
 (* ev->expr[j]) (&i, &t, ev);
 if (i == -123456 && t == - TEND_EVENT) {

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
    ev->i = -1; ev->t = - TEND_EVENT;
    if (ev->expr[0]) {
      (* ev->expr[0]) (&ev->i, &ev->t, ev);
      if (ev->i == END_EVENT || ev->t == TEND_EVENT) {
 ev->i = END_EVENT; ev->t = - TEND_EVENT;
      }
    }
    else if (ev->expr[2]) {
      (* ev->expr[2]) (&ev->i, &ev->t, ev);
      if (ev->i != -1)
 ev->i = 0;
      if (ev->t != - TEND_EVENT)
 ev->t = 0;
    }
  }
}

enum { event_done, event_alive, event_stop };

static int event_finished (Event * ev)
{
  ev->i = -1; ev->t = - TEND_EVENT;
  return event_done;
}

void event_register (Event event) {
  if (!(Events)) qassert ("/home/lisergey/basilisk/src/grid/events.h", 88, "Events");
  if (!(!event.last)) qassert ("/home/lisergey/basilisk/src/grid/events.h", 89, "!event.last");
  int n = 0, parent = -1;
  for (Event * ev = Events; !ev->last; ev++) {
    if (!strcmp (event.name, ev->name)) {
      if (!(parent < 0)) qassert ("/home/lisergey/basilisk/src/grid/events.h", 93, "parent < 0");
      parent = n;
    }
    n++;
  }
  if (parent < 0) {
    Events = (Event *) prealloc (Events, (n + 2)*sizeof(Event),__func__,__FILE__,__LINE__);
    Events[n] = event;
    Events[n].next = NULL;
    Events[n + 1].last = true;
    init_event (&Events[n]);
  }
  else {
    Event * ev = ((Event *) pcalloc (1, sizeof(Event),__func__,__FILE__,__LINE__));
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
#line 136 "/home/lisergey/basilisk/src/grid/events.h"
static bool overload_event() { return true; }

static int event_do (Event * ev, bool action)
{
  if ((iter > ev->i && t > ev->t) || !event_cond (ev, iter, t))
    return event_finished (ev);
  if (!overload_event() || iter == ev->i || fabs (t - ev->t) <= TEPS*t) {
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
    if (ev->i == END_EVENT && action)
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
  inext = END_EVENT; tnext = 1e30;
  for (Event * ev = Events; !ev->last && !cond; ev++)
    if (ev->i != END_EVENT &&
 (ev->expr[1] || (ev->expr[0] && !ev->expr[1] && !ev->expr[2]) || ev->arrayi || ev->arrayt))
      cond = 1;
  for (Event * ev = Events; !ev->last; ev++) {
    int status = event_do (ev, action);
    if (status == event_stop) {
      end_event_do (action);
      return 0;
    }
    if (status == event_alive && ev->i != END_EVENT &&
 (ev->expr[1] || (ev->expr[0] && !ev->expr[1] && !ev->expr[2]) || ev->arrayi || ev->arrayt))
      cond1 = 1;
    if (ev->t > t && ev->t < tnext)
      tnext = ev->t;
    if (ev->i > iter && ev->i < inext)
      inext = ev->i;
  }
  if (overload_event() && (!cond || cond1) && (tnext != 1e30 || inext != END_EVENT)) {
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
  if (tnext != 1e30 && tnext > t) {
    unsigned int n = (tnext - t)/dt;
    if (!(n < INT_MAX)) qassert ("/home/lisergey/basilisk/src/grid/events.h", 252, "n < INT_MAX");
    if (n == 0)
      dt = tnext - t;
    else {
      double dt1 = (tnext - t)/n;
      if (dt1 > dt*(1. + TEPS))
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
#line 2 "/home/lisergey/basilisk/src/grid/cartesian-common.h"

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
  double Delta_z = Delta;\
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
  double z = (kg/2. + (point.k - 2) + 0.5)*Delta + Z0;\
\
\
\
  NOT_UNUSED(z);\
\
  NOT_UNUSED(Delta);\
  NOT_UNUSED(Delta_x);\
\
  NOT_UNUSED(Delta_y);\
\
\
  NOT_UNUSED(Delta_z);\
\
\
  ;\

#line 44


#line 1 "grid/fpe.h"
#line 1 "/home/lisergey/basilisk/src/grid/fpe.h"


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
#line 47 "/home/lisergey/basilisk/src/grid/cartesian-common.h"

#define end_foreach_face()

#define foreach_point(...)\
{\
  int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);\
  coord _p = { __VA_ARGS__ };\
  Point point = locate (_p.x, _p.y, _p.z);\
  if (point.level >= 0) {\
    POINT_VARIABLES\

#line 57

#define end_foreach_point() }}

#define foreach_region(p, box, n)\
  OMP_PARALLEL() { NOT_UNUSED (p);\
    coord p = {0, 0, box[0].z};\
    OMP(omp for schedule(static))\
      for (int _i = 0; _i < (int) n.x; _i++) {\
 p.x = box[0].x + (box[1].x - box[0].x)/n.x*(_i + 0.5);\
 for (int _j = 0; _j < (int) n.y; _j++) {\
   p.y = box[0].y + (box[1].y - box[0].y)/n.y*(_j + 0.5);\
   Point point = locate (p.x, p.y, p.z);\
   if (point.level >= 0) {\
     int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);\
     POINT_VARIABLES\

#line 72

#define end_foreach_region() }}}}





static void init_block_scalar (scalar sb, const char * name, const char * ext,
          int n, int block)
{
  char bname[strlen(name) + strlen(ext) + 10];
  if (n == 0) {
    sprintf (bname, "%s%s", name, ext);
    _attribute[sb.i].block = block;
    baseblock = list_append (baseblock, sb);
  }
  else {
    sprintf (bname, "%s%d%s", name, n, ext);
    _attribute[sb.i].block = - n;
  }
  _attribute[sb.i].name = pstrdup (bname,__func__,__FILE__,__LINE__);
  all = list_append (all, sb);
}

#define interpreter_set_int(...)
#define interpreter_reset_scalar(...)

scalar alloc_block_scalar (const char * name, const char * ext, int block)
{
  interpreter_set_int (&block);
  int nvar = datasize/sizeof(double);

  scalar s = {0};
  while (s.i < nvar) {
    int n = 0;
    scalar sb = s;
    while (sb.i < nvar && n < block && _attribute[sb.i].freed)
      n++, sb.i++;
    if (n >= block) {
      memset (&_attribute[s.i], 0, block*sizeof (_Attributes));
      for (sb.i = s.i, n = 0; n < block; n++, sb.i++) {
 init_block_scalar (sb, name, ext, n, block);
 interpreter_reset_scalar (sb);
      }
      trash (((scalar []){s, {-1}}));
      return s;
    }
    s.i = sb.i + 1;
  }


  s = (scalar){nvar};
  if (!(nvar + block <= _NVARMAX)) qassert ("/home/lisergey/basilisk/src/grid/cartesian-common.h", 124, "nvar + block <= _NVARMAX");

  if (_attribute == NULL)
    _attribute = (_Attributes *) pcalloc (nvar + block + 1, sizeof (_Attributes),__func__,__FILE__,__LINE__);
  else
    _attribute = (_Attributes *)
      prealloc (_attribute, (nvar + block + 1)*sizeof (_Attributes),__func__,__FILE__,__LINE__);
  memset (&_attribute[nvar], 0, block*sizeof (_Attributes));
  for (int n = 0; n < block; n++, nvar++) {
    scalar sb = (scalar){nvar};
    init_block_scalar (sb, name, ext, n, block);
  }

  realloc_scalar (block*sizeof(double));
  trash (((scalar []){s, {-1}}));
  return s;
}

scalar new_block_scalar (const char * name, const char * ext, int block)
{
  scalar s = alloc_block_scalar (name, ext, block), sb;
  int n = 0;
  for (sb.i = s.i, n = 0; n < block; n++, sb.i++)
    init_scalar (sb, NULL);
  return s;
}

scalar new_scalar (const char * name)
{
  return init_scalar (alloc_block_scalar (name, "", 1), NULL);
}

scalar new_vertex_scalar (const char * name)
{
  return init_vertex_scalar (alloc_block_scalar (name, "", 1), NULL);
}

static vector alloc_block_vector (const char * name, int block)
{
  vector v;
  struct { char * x, * y, * z; } ext = {".x", ".y", ".z"};
  
    v.x = alloc_block_scalar (name, ext.x, block);
    
#line 166
v.y = alloc_block_scalar (name, ext.y, block);
    
#line 166
v.z = alloc_block_scalar (name, ext.z, block);
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
      
#line 190
vb.y.i = v.y.i + i;
      
#line 190
vb.z.i = v.z.i + i;
    init_vector (vb, NULL);
    
      _attribute[vb.x.i].block = - i;
      
#line 193
_attribute[vb.y.i].block = - i;
      
#line 193
_attribute[vb.z.i].block = - i;
  }
  
    _attribute[v.x.i].block = block;
    
#line 196
_attribute[v.y.i].block = block;
    
#line 196
_attribute[v.z.i].block = block;
  return v;
}

vector new_block_face_vector (const char * name, int block)
{
  vector v = alloc_block_vector (name, block);
  for (int i = 0; i < block; i++) {
    vector vb;
    
      vb.x.i = v.x.i + i;
      
#line 206
vb.y.i = v.y.i + i;
      
#line 206
vb.z.i = v.z.i + i;
    init_face_vector (vb, NULL);
    
      _attribute[vb.x.i].block = - i;
      
#line 209
_attribute[vb.y.i].block = - i;
      
#line 209
_attribute[vb.z.i].block = - i;
  }
  
    _attribute[v.x.i].block = block;
    
#line 212
_attribute[v.y.i].block = block;
    
#line 212
_attribute[v.z.i].block = block;
  return v;
}

tensor new_tensor (const char * name)
{
  char cname[strlen(name) + 3];
  struct { char * x, * y, * z; } ext = {"%s.x", "%s.y", "%s.z"};
  tensor t;
   {
    sprintf (cname, ext.x, name);
    t.x = alloc_block_vector (cname, 1);
  } 
#line 221
{
    sprintf (cname, ext.y, name);
    t.y = alloc_block_vector (cname, 1);
  } 
#line 221
{
    sprintf (cname, ext.z, name);
    t.z = alloc_block_vector (cname, 1);
  }
  init_tensor (t, NULL);
  return t;
}

tensor new_symmetric_tensor (const char * name)
{
  struct { char * x, * y, * z; } ext = {".x.x", ".y.y", ".z.z"};
  tensor t;
  
    t.x.x = alloc_block_scalar (name, ext.x, 1);
    
#line 234
t.y.y = alloc_block_scalar (name, ext.y, 1);
    
#line 234
t.z.z = alloc_block_scalar (name, ext.z, 1);

    t.x.y = alloc_block_scalar (name, ".x.y", 1);
    t.y.x = t.x.y;


    t.x.z = alloc_block_scalar (name, ".x.z", 1);
    t.z.x = t.x.z;
    t.y.z = alloc_block_scalar (name, ".y.z", 1);
    t.z.y = t.y.z;




  init_tensor (t, NULL);
  return t;
}

static int nconst = 0;

void init_const_scalar (scalar s, const char * name, double val)
{
  if (s.i - _NVARMAX >= nconst) {
    nconst = s.i - _NVARMAX + 1;
    _constant = (double *) prealloc (_constant, (nconst)*sizeof(double),__func__,__FILE__,__LINE__);
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
    
#line 273
init_const_scalar (v.y, name, *val++);
    
#line 273
init_const_scalar (v.z, name, *val++);
}

vector new_const_vector (const char * name, int i, double * val)
{
  vector v;
  
    v.x.i = _NVARMAX + i++;
    
#line 280
v.y.i = _NVARMAX + i++;
    
#line 280
v.z.i = _NVARMAX + i++;
  init_const_vector (v, name, val);
  return v;
}

static void cartesian_scalar_clone (scalar clone, scalar src)
{
  char * cname = _attribute[clone.i].name;
  double (** boundary) (Point, Point, scalar, void *) = _attribute[clone.i].boundary;
  double (** boundary_homogeneous) (Point, Point, scalar, void *) =
    _attribute[clone.i].boundary_homogeneous;
  if (!(_attribute[src.i].block > 0 && _attribute[clone.i].block == _attribute[src.i].block)) qassert ("/home/lisergey/basilisk/src/grid/cartesian-common.h", 291, "src.block > 0 && clone.block == src.block");
  pfree (_attribute[clone.i].depends,__func__,__FILE__,__LINE__);
  _attribute[clone.i] = _attribute[src.i];
  _attribute[clone.i].name = cname;
  _attribute[clone.i].boundary = boundary;
  _attribute[clone.i].boundary_homogeneous = boundary_homogeneous;
  for (int i = 0; i < nboundary; i++) {
    _attribute[clone.i].boundary[i] = _attribute[src.i].boundary[i];
    _attribute[clone.i].boundary_homogeneous[i] = _attribute[src.i].boundary_homogeneous[i];
  }
  _attribute[clone.i].depends = list_copy (_attribute[src.i].depends);
}

scalar * list_clone (scalar * l)
{
  scalar * list = NULL;
  int nvar = datasize/sizeof(double), map[nvar];
  for (int i = 0; i < nvar; i++)
    map[i] = -1;
  {scalar*_i=(scalar*)( l);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {
    scalar c = _attribute[s.i].block > 1 ? new_block_scalar("c", "", _attribute[s.i].block) : new_scalar("c");
    scalar_clone (c, s);
    map[s.i] = c.i;
    list = list_append (list, c);
  }}}
  {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    {
      if (_attribute[s.i].v.x.i >= 0 && map[_attribute[s.i].v.x.i] >= 0)
 _attribute[s.i].v.x.i = map[_attribute[s.i].v.x.i];
      
#line 318
if (_attribute[s.i].v.y.i >= 0 && map[_attribute[s.i].v.y.i] >= 0)
 _attribute[s.i].v.y.i = map[_attribute[s.i].v.y.i];
      
#line 318
if (_attribute[s.i].v.z.i >= 0 && map[_attribute[s.i].v.z.i] >= 0)
 _attribute[s.i].v.z.i = map[_attribute[s.i].v.z.i];}}}
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
      pfree (_attribute[fb.i].name,__func__,__FILE__,__LINE__); _attribute[fb.i].name = NULL;
      pfree (_attribute[fb.i].boundary,__func__,__FILE__,__LINE__); _attribute[fb.i].boundary = NULL;
      pfree (_attribute[fb.i].boundary_homogeneous,__func__,__FILE__,__LINE__); _attribute[fb.i].boundary_homogeneous = NULL;
      pfree (_attribute[fb.i].depends,__func__,__FILE__,__LINE__); _attribute[fb.i].depends = NULL;
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
      scalar * s;
      for (s = all; s->i >= 0 && s->i != f.i; s++);
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
  if (!(_val_higher_dimension == 0.)) qassert ("/home/lisergey/basilisk/src/grid/cartesian-common.h", 369, "_val_higher_dimension == 0.");

  if (free_solver_funcs) {
    free_solver_func * a = (free_solver_func *) free_solver_funcs->p;
    for (int i = 0; i < free_solver_funcs->len/sizeof(free_solver_func); i++)
      a[i] ();
    array_free (free_solver_funcs);
  }

  delete (all);
  pfree (all,__func__,__FILE__,__LINE__); all = NULL;
  pfree (baseblock,__func__,__FILE__,__LINE__); baseblock = NULL;
  for (Event * ev = Events; !ev->last; ev++) {
    Event * e = ev->next;
    while (e) {
      Event * next = e->next;
      pfree (e,__func__,__FILE__,__LINE__);
      e = next;
    }
  }

  pfree (Events,__func__,__FILE__,__LINE__); Events = NULL;
  pfree (_attribute,__func__,__FILE__,__LINE__); _attribute = NULL;
  pfree (_constant,__func__,__FILE__,__LINE__); _constant = NULL;
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
      
#line 421
list1.y = list_append (list1.y, v.y);
      
#line 421
list1.z = list_append (list1.z, v.z);}}}
  boundary_face (list1);
  
    pfree (list1.x,__func__,__FILE__,__LINE__);
    
#line 424
pfree (list1.y,__func__,__FILE__,__LINE__);
    
#line 424
pfree (list1.z,__func__,__FILE__,__LINE__);
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
{tracing("boundary_internal","/home/lisergey/basilisk/src/grid/cartesian-common.h",440);
  if (list == NULL)
    {end_tracing("boundary_internal","/home/lisergey/basilisk/src/grid/cartesian-common.h",443);return;}
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
     
#line 452
if (_attribute[s.i].v.y.i == s.i)
       listf.y = list_add (listf.y, s), flux = true;
     
#line 452
if (_attribute[s.i].v.z.i == s.i)
       listf.z = list_add (listf.z, s), flux = true;}
 if (!is_constant(cm) && _attribute[cm.i].dirty)
   listc = list_add_depends (listc, cm);
 if (_attribute[s.i].face != 2)
   listc = list_add_depends (listc, s);
      }




    }}}
  if (flux) {
    boundary_face (listf);
    
      pfree (listf.x,__func__,__FILE__,__LINE__);
      
#line 467
pfree (listf.y,__func__,__FILE__,__LINE__);
      
#line 467
pfree (listf.z,__func__,__FILE__,__LINE__);
  }
  if (listc) {
    boundary_level (listc, -1);
    {scalar*_i=(scalar*)( listc);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      _attribute[s.i].dirty = false;}}
    pfree (listc,__func__,__FILE__,__LINE__);
  }
end_tracing("boundary_internal","/home/lisergey/basilisk/src/grid/cartesian-common.h",475);}

void cartesian_boundary_level (scalar * list, int l)
{
  { Boundary ** _i = boundaries, * _b; while (_i && (_b = *_i++)) if (_b->level) _b->level (_b, list, l); };
}

void cartesian_boundary_face (vectorl list)
{
  
    {scalar*_i=(scalar*)( list.x);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      _attribute[s.i].dirty = 2;}}
    
#line 485
{scalar*_i=(scalar*)( list.y);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      _attribute[s.i].dirty = 2;}}
    
#line 485
{scalar*_i=(scalar*)( list.z);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      _attribute[s.i].dirty = 2;}}
}

static double symmetry (Point point, Point neighbor, scalar s, void * data)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;
  return val(s,0,0,0);
}

static double antisymmetry (Point point, Point neighbor, scalar s, void * data)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;
  return -val(s,0,0,0);
}

double (* default_scalar_bc[]) (Point, Point, scalar, void *) = {
  symmetry, symmetry, symmetry, symmetry, symmetry, symmetry
};

scalar cartesian_init_scalar (scalar s, const char * name)
{

  char * pname;
  if (name) {
    pfree (_attribute[s.i].name,__func__,__FILE__,__LINE__);
    pname = pstrdup (name,__func__,__FILE__,__LINE__);
  }
  else
    pname = _attribute[s.i].name;
  int block = _attribute[s.i].block;
  double (** boundary) (Point, Point, scalar, void *) = _attribute[s.i].boundary;
  double (** boundary_homogeneous) (Point, Point, scalar, void *) =
    _attribute[s.i].boundary_homogeneous;
  _attribute[s.i].name = pname;
  if (block < 0)
    _attribute[s.i].block = block;
  else
    _attribute[s.i].block = block > 0 ? block : 1;

  _attribute[s.i].boundary = boundary ? boundary :
    (double (**)(Point, Point, scalar, void *))
    pmalloc (nboundary*sizeof (void (*)()),__func__,__FILE__,__LINE__);
  _attribute[s.i].boundary_homogeneous = boundary_homogeneous ? boundary_homogeneous :
    (double (**)(Point, Point, scalar, void *))
    pmalloc (nboundary*sizeof (void (*)()),__func__,__FILE__,__LINE__);
  for (int b = 0; b < nboundary; b++)
    _attribute[s.i].boundary[b] = _attribute[s.i].boundary_homogeneous[b] =
      b < 2*3 ? default_scalar_bc[b] : symmetry;
  _attribute[s.i].gradient = NULL;
   {
    _attribute[s.i].d.x = 0;
    _attribute[s.i].v.x.i = -1;
  } 
#line 533
{
    _attribute[s.i].d.y = 0;
    _attribute[s.i].v.y.i = -1;
  } 
#line 533
{
    _attribute[s.i].d.z = 0;
    _attribute[s.i].v.z.i = -1;
  }
  _attribute[s.i].face = false;
  return s;
}

scalar cartesian_init_vertex_scalar (scalar s, const char * name)
{
  cartesian_init_scalar (s, name);
  
    _attribute[s.i].d.x = -1;
    
#line 545
_attribute[s.i].d.y = -1;
    
#line 545
_attribute[s.i].d.z = -1;
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
      cartesian_init_scalar (v.x, cname);
    }
    else
      cartesian_init_scalar (v.x, NULL);
    _attribute[v.x.i].v = v;
  } 
#line 560
{
    if (name) {
      char cname[strlen(name) + 3];
      sprintf (cname, "%s%s", name, ext.y);
      cartesian_init_scalar (v.y, cname);
    }
    else
      cartesian_init_scalar (v.y, NULL);
    _attribute[v.y.i].v = v;
  } 
#line 560
{
    if (name) {
      char cname[strlen(name) + 3];
      sprintf (cname, "%s%s", name, ext.z);
      cartesian_init_scalar (v.z, cname);
    }
    else
      cartesian_init_scalar (v.z, NULL);
    _attribute[v.z.i].v = v;
  }

  for (int d = 0; d < nboundary; d++)
    _attribute[v.x.i].boundary[d] = _attribute[v.x.i].boundary_homogeneous[d] =
      d < 2*3 ? default_vector_bc[d] : antisymmetry;
  return v;
}

vector cartesian_init_face_vector (vector v, const char * name)
{
  v = cartesian_init_vector (v, name);
   {
    _attribute[v.x.i].d.x = -1;
    _attribute[v.x.i].face = true;
  } 
#line 580
{
    _attribute[v.y.i].d.y = -1;
    _attribute[v.y.i].face = true;
  } 
#line 580
{
    _attribute[v.z.i].d.z = -1;
    _attribute[v.z.i].face = true;
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
      cartesian_init_vector (t.x, cname);
    }
    else
      cartesian_init_vector (t.x, NULL);
  } 
#line 592
{
    if (name) {
      char cname[strlen(name) + 3];
      sprintf (cname, "%s%s", name, ext.y);
      cartesian_init_vector (t.y, cname);
    }
    else
      cartesian_init_vector (t.y, NULL);
  } 
#line 592
{
    if (name) {
      char cname[strlen(name) + 3];
      sprintf (cname, "%s%s", name, ext.z);
      cartesian_init_vector (t.z, cname);
    }
    else
      cartesian_init_vector (t.z, NULL);
  }
#line 616 "/home/lisergey/basilisk/src/grid/cartesian-common.h"
    if (!(false)) qassert ("/home/lisergey/basilisk/src/grid/cartesian-common.h", 616, "false");

  return t;
}

void output_cells (FILE * fp, coord c, double size)
{
  {foreach() {
    bool inside = true;
    coord o = {x,y,z};
    
      if (inside && size > 0. &&
   (o.x > c.x + size || o.x < c.x - size))
 inside = false;
      
#line 627
if (inside && size > 0. &&
   (o.y > c.y + size || o.y < c.y - size))
 inside = false;
      
#line 627
if (inside && size > 0. &&
   (o.z > c.z + size || o.z < c.z - size))
 inside = false;
    if (inside) {
      Delta /= 2.;
#line 642 "/home/lisergey/basilisk/src/grid/cartesian-common.h"
      for (int i = -1; i <= 1; i += 2) {
 fprintf (fp, "%g %g %g\n%g %g %g\n%g %g %g\n%g %g %g\n%g %g %g\n\n",
   x - Delta, y - Delta, z + i*Delta,
   x - Delta, y + Delta, z + i*Delta,
   x + Delta, y + Delta, z + i*Delta,
   x + Delta, y - Delta, z + i*Delta,
   x - Delta, y - Delta, z + i*Delta);
 for (int j = -1; j <= 1; j += 2)
   fprintf (fp, "%g %g %g\n%g %g %g\n\n",
     x + i*Delta, y + j*Delta, z - Delta,
     x + i*Delta, y + j*Delta, z + Delta);
      }

    }
  }end_foreach();}
  fflush (fp);
}


static void output_cells_internal (FILE * fp)
{
  output_cells (fp
#line 621
,( coord) {0}, 0.
#line 663
);
}


static char * replace_ (const char * vname)
{
  char * name = pstrdup (vname,__func__,__FILE__,__LINE__), * c = name;
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







    "  splot '%s' w l lc 0, "
    "'%s' u 1+4*v:2+4*v:3+4*v:4+4*v w labels tc lt 1"
           " title columnhead(4+4*v)",

    vname, cells, stencil);
  pfree (vname,__func__,__FILE__,__LINE__);
}

void cartesian_debug (Point point)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;
  char name[80] = "cells";
  if (pid() > 0)
    sprintf (name, "cells-%d", pid());
  FILE * fp = fopen (name, "w");
  output_cells (fp, (coord){x,y,z}, 4.*Delta);
  fclose (fp);

  char stencil[80] = "stencil";
  if (pid() > 0)
    sprintf (stencil, "stencil-%d", pid());
  fp = fopen (stencil, "w");
  {scalar*_i=(scalar*)( all);if(_i)for(scalar v=*_i;(&v)->i>=0;v=*++_i){





    fprintf (fp, "x y z %s ", _attribute[v.i].name);}}

  fputc ('\n', fp);
#line 748 "/home/lisergey/basilisk/src/grid/cartesian-common.h"
    for (int k = -2; k <= 2; k++)
      for (int l = -2; l <= 2; l++)
 for (int m = -2; m <= 2; m++) {
   {scalar*_i=(scalar*)( all);if(_i)for(scalar v=*_i;(&v)->i>=0;v=*++_i){ {
     fprintf (fp, "%g %g %g ",
       x + k*Delta + _attribute[v.i].d.x*Delta/2.,
       y + l*Delta + _attribute[v.i].d.y*Delta/2.,
       z + m*Delta + _attribute[v.i].d.z*Delta/2.);
     if (allocated(k,l,m))
       fprintf (fp, "%g ", val(v,k,l,m));
     else
       fputs ("n/a ", fp);
   }}}
   fputc ('\n', fp);
 }

  fclose (fp);

  fp = fopen ("debug.plot", "w");
  fprintf (fp,
    "set term x11\n"
    "set size ratio -1\n"
    "set key outside\n");
  {scalar*_i=(scalar*)( all);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {
    char * name = replace_ (_attribute[s.i].name);
    fprintf (fp, "%s = %d\n", name, s.i);
    pfree (name,__func__,__FILE__,__LINE__);
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
  init_face_vector = cartesian_init_face_vector;
  init_tensor = cartesian_init_tensor;
  boundary_level = cartesian_boundary_level;
  boundary_face = cartesian_boundary_face;
  scalar_clone = cartesian_scalar_clone;
  debug = cartesian_debug;
}

tensor init_symmetric_tensor (tensor t, const char * name)
{
  return init_tensor (t, name);
}

static double interpolate_linear (Point point, scalar v,
      double xp, double yp, double zp)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;
#line 823 "/home/lisergey/basilisk/src/grid/cartesian-common.h"
  x = (xp - x)/Delta - _attribute[v.i].d.x/2.;
  y = (yp - y)/Delta - _attribute[v.i].d.y/2.;
  z = (zp - z)/Delta - _attribute[v.i].d.z/2.;
  int i = sign(x), j = sign(y), k = sign(z);
  x = fabs(x); y = fabs(y); z = fabs(z);

  return (((val(v,0,0,0)*(1. - x) + val(v,i,0,0)*x)*(1. - y) +
    (val(v,0,j,0)*(1. - x) + val(v,i,j,0)*x)*y)*(1. - z) +
   ((val(v,0,0,k)*(1. - x) + val(v,i,0,k)*x)*(1. - y) +
    (val(v,0,j,k)*(1. - x) + val(v,i,j,k)*x)*y)*z);

}


#line 805
static void _stencil_interpolate_linear (Point point, scalar v,
_stencil_undefined * xp,_stencil_undefined * yp,_stencil_undefined * zp)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;         
#line 823 "/home/lisergey/basilisk/src/grid/cartesian-common.h"
        
        
        
  
          

_stencil_val(v,0,0,0);_stencil_val(v, o_stencil,0,0);
_stencil_val(v,0,o_stencil,0); _stencil_val(v,o_stencil,o_stencil,0);
_stencil_val(v,0,0,o_stencil); _stencil_val(v,o_stencil,0,o_stencil);
_stencil_val(v,0,o_stencil,o_stencil); _stencil_val(v,o_stencil,o_stencil,o_stencil);

  
#line 829
return           
              
    
    ;

}

     
double interpolate (scalar v, double xp, double yp, double zp,
      bool linear)
{tracing("interpolate","/home/lisergey/basilisk/src/grid/cartesian-common.h",837);
  double val = 1e30;
  foreach_point_stencil (1,{(NonLocal[]){{"v","scalar",(void *)&v,NULL,0},{"linear","bool",(void *)&linear,NULL,0},{"val","double",(void *)&val,NULL,0,'m'},{"zp","double",(void *)&zp,NULL,0},{"yp","double",(void *)&yp,NULL,0},{"xp","double",(void *)&xp,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"\n\n\n// #line 805 \"/home/lisergey/basilisk/src/grid/cartesian-common.h\"\nstatic real interpolate_linear (Point point, scalar v,\n      real xp, real yp, real zp)\n{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;\n// # 823 \"/home/lisergey/basilisk/src/grid/cartesian-common.h\"\n  x = (xp - x)/Delta - v.d.x/2.;\n  y = (yp - y)/Delta - v.d.y/2.;\n  z = (zp - z)/Delta - v.d.z/2.;\n  int i = sign(x), j = sign(y), k = sign(z);\n  x = fabs(x); y = fabs(y); z = fabs(z);\n\n  return (((val(v,0,0,0)*(1. - x) + val(v,i,0,0)*x)*(1. - y) +\n    (val(v,0,j,0)*(1. - x) + val(v,i,j,0)*x)*y)*(1. - z) +\n   ((val(v,0,0,k)*(1. - x) + val(v,i,0,k)*x)*(1. - y) +\n    (val(v,0,j,k)*(1. - x) + val(v,i,j,k)*x)*y)*z);\n\n}","\n    \n// #line 842 \"/home/lisergey/basilisk/src/grid/cartesian-common.h\"\nval = linear ? interpolate_linear (point, v, xp, yp, zp) : val(v,0,0,0);"})
    { _stencil_interpolate_linear (point, v, NULL, NULL, NULL); _stencil_val(v,0,0,0);    }end_foreach_point_stencil();
  
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel  reduction (min:val)){
#line 841
foreach_point (xp, yp, zp)
    val = linear ? interpolate_linear (point, v, xp, yp, zp) : val(v,0,0,0);end_foreach_point();mpi_all_reduce_array(&val,double,MPI_MIN,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
  
#line 843
{end_tracing("interpolate","/home/lisergey/basilisk/src/grid/cartesian-common.h",843);return val;}
end_tracing("interpolate","/home/lisergey/basilisk/src/grid/cartesian-common.h",844);}

     
void interpolate_array (scalar * list, coord * a, int n, double * v,
   bool linear)
{tracing("interpolate_array","/home/lisergey/basilisk/src/grid/cartesian-common.h",847);
  int len = 0;
  {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    len++;}}
  for (int i = 0; i < n; i++) {
    double * w = v;
    {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      *(w++) = 1e30;}}
    foreach_point_stencil (1,{(NonLocal[]){{"linear","bool",(void *)&linear,NULL,0},{"v","double",(void *)v,NULL,1},{"list","scalar",(void *)list,NULL,1},{"len","int",(void *)&len,NULL,0},{"i","int",(void *)&i,NULL,0},{"a","coord",(void *)a,NULL,1},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"\n\n\n// #line 805 \"/home/lisergey/basilisk/src/grid/cartesian-common.h\"\nstatic real interpolate_linear (Point point, scalar v,\n      real xp, real yp, real zp)\n{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;\n// # 823 \"/home/lisergey/basilisk/src/grid/cartesian-common.h\"\n  x = (xp - x)/Delta - v.d.x/2.;\n  y = (yp - y)/Delta - v.d.y/2.;\n  z = (zp - z)/Delta - v.d.z/2.;\n  int i = sign(x), j = sign(y), k = sign(z);\n  x = fabs(x); y = fabs(y); z = fabs(z);\n\n  return (((val(v,0,0,0)*(1. - x) + val(v,i,0,0)*x)*(1. - y) +\n    (val(v,0,j,0)*(1. - x) + val(v,i,j,0)*x)*y)*(1. - z) +\n   ((val(v,0,0,k)*(1. - x) + val(v,i,0,k)*x)*(1. - y) +\n    (val(v,0,j,k)*(1. - x) + val(v,i,j,k)*x)*y)*z);\n\n}"," \n// #line 857 \"/home/lisergey/basilisk/src/grid/cartesian-common.h\"\n{\n      int j = 0;\n      {forin (scalar, s , list)\n v[j++] = !linear ? val(s,0,0,0) : interpolate_linear (point, s, a[i].x, a[i].y, a[i].z); endforin()}\n    }"}) {   
      
      {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
 { _stencil_val(s,0,0,0); _stencil_interpolate_linear (point, s, NULL, NULL, NULL);    }}}
    }end_foreach_point_stencil();
    
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel  reduction(min:v[:len])){
#line 857
foreach_point (a[i].x, a[i].y, a[i].z) {
      int j = 0;
      {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
 v[j++] = !linear ? val(s,0,0,0) : interpolate_linear (point, s, a[i].x, a[i].y, a[i].z);}}
    }end_foreach_point();mpi_all_reduce_array(v,double,MPI_MIN,len);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
    
#line 862
v = w;
  }
end_tracing("interpolate_array","/home/lisergey/basilisk/src/grid/cartesian-common.h",864);}



typedef int bid;

bid new_bid()
{
  int b = nboundary++;
  {scalar*_i=(scalar*)( all);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {
    _attribute[s.i].boundary = (double (**)(Point, Point, scalar, void *))
      prealloc (_attribute[s.i].boundary, nboundary*sizeof (void (*)()),__func__,__FILE__,__LINE__);
    _attribute[s.i].boundary_homogeneous = (double (**)(Point, Point, scalar, void *))
      prealloc (_attribute[s.i].boundary_homogeneous, nboundary*sizeof (void (*)()),__func__,__FILE__,__LINE__);
  }}}
  {scalar*_i=(scalar*)( all);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {
    if (_attribute[s.i].v.x.i < 0)
      _attribute[s.i].boundary[b] = _attribute[s.i].boundary_homogeneous[b] = symmetry;
    else if (_attribute[s.i].v.x.i == s.i) {
      vector v = _attribute[s.i].v;
      
 _attribute[v.y.i].boundary[b] = _attribute[v.y.i].boundary_homogeneous[b] = symmetry;
 
#line 885
_attribute[v.z.i].boundary[b] = _attribute[v.z.i].boundary_homogeneous[b] = symmetry;
 
#line 885
_attribute[v.x.i].boundary[b] = _attribute[v.x.i].boundary_homogeneous[b] = symmetry;
      _attribute[v.x.i].boundary[b] = _attribute[v.x.i].boundary_homogeneous[b] =
 _attribute[v.x.i].face ? NULL : antisymmetry;
    }
  }}}
  return b;
}



static double periodic_bc (Point point, Point neighbor, scalar s, void * data)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;
  return 1e30;
}

static void periodic_boundary (int d)
{

  {scalar*_i=(scalar*)( all);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    if (is_vertex_scalar (s))
      _attribute[s.i].boundary[d] = _attribute[s.i].boundary_homogeneous[d] = NULL;
    else
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





    if (!(dir <= back)) qassert ("/home/lisergey/basilisk/src/grid/cartesian-common.h", 926, "dir <= back");


  int c = dir/2;
  periodic_boundary (2*c);
  periodic_boundary (2*c + 1);
  (&Period.x)[c] = true;
}


double getvalue (Point point, scalar s, int i, int j, int k)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;
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
  for (int d = 1; d < 3; d++)
    fprintf (qstderr(), ",%d", index[d]);
  fputs ("]", qstderr());
}

void stencil_val (Point p, scalar s, int i, int j, int k,
    const char * file, int line, bool overflow)
{
  if (is_constant(s) || s.i < 0)
    return;
  int index[] = {i, j, k};
  for (int d = 0; d < 3; d++)
    index[d] += (&p.i)[d];
  bool central = true;
  for (int d = 0; d < 3; d++) {
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
#line 986
{
      if ((!_attribute[s.i].face || _attribute[s.i].v.y.i != s.i) && abs(index[d]) > _attribute[s.i].width)
 _attribute[s.i].width = abs(index[d]);
      d++;
    } 
#line 986
{
      if ((!_attribute[s.i].face || _attribute[s.i].v.z.i != s.i) && abs(index[d]) > _attribute[s.i].width)
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
  for (int d = 0; d < 3; d++)
    index[d] += (&p.i)[d];
  for (int d = 0; d < 3; d++)
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




#define dimensional(...)

#define show_dimension_internal(...)
#define display_value(...)
#define interpreter_verbosity(...)
#line 4 "/home/lisergey/basilisk/src/grid/multigrid-common.h"

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
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;
  double sum = 0.;
  {foreach_child()
    sum += val(s,0,0,0);end_foreach_child()}
  val(s,0,0,0) = sum/(1 << 3);
}

static inline void restriction_volume_average (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;

#line 35
if(!is_constant(cm)){{
  double sum = 0.;
  {foreach_child()
    sum += val(cm,0,0,0)*val(s,0,0,0);end_foreach_child()}
  val(s,0,0,0) = sum/(1 << 3)/(val(cm,0,0,0) + 1e-30);
}}else {double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);

#line 35
{
  double sum = 0.;
  {foreach_child()
    sum += _const_cm*val(s,0,0,0);end_foreach_child()}
  val(s,0,0,0) = sum/(1 << 3)/(_const_cm + 1e-30);
}}

#line 40
}

static inline void face_average (Point point, vector v)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;
   {







      val(v.x,0,0,0) = (fine(v.x,0,0,0) + fine(v.x,0,1,0) +
        fine(v.x,0,0,1) + fine(v.x,0,1,1))/4.;
      val(v.x,1,0,0) = (fine(v.x,2,0,0) + fine(v.x,2,1,0) +
  fine(v.x,2,0,1) + fine(v.x,2,1,1))/4.;

  } 
#line 44
{







      val(v.y,0,0,0) = (fine(v.y,0,0,0) + fine(v.y,0,0,1) +
        fine(v.y,1,0,0) + fine(v.y,1,0,1))/4.;
      val(v.y,0,1,0) = (fine(v.y,0,2,0) + fine(v.y,0,2,1) +
  fine(v.y,1,2,0) + fine(v.y,1,2,1))/4.;

  } 
#line 44
{







      val(v.z,0,0,0) = (fine(v.z,0,0,0) + fine(v.z,1,0,0) +
        fine(v.z,0,1,0) + fine(v.z,1,1,0))/4.;
      val(v.z,0,0,1) = (fine(v.z,0,0,2) + fine(v.z,1,0,2) +
  fine(v.z,0,1,2) + fine(v.z,1,1,2))/4.;

  }
}

static inline void restriction_face (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;
  face_average (point, _attribute[s.i].v);
}

static inline void restriction_vertex (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;
  for (int i = 0; i <= 1; i++) {
    val(s,i,0,0) = fine(s,2*i,0,0);

    val(s,i,1,0) = fine(s,2*i,2,0);


    for (int j = 0; j <= 1; j++)
      val(s,i,j,1) = fine(s,2*i,2*j,2);

  }
}

static inline void no_restriction (Point point, scalar s) {int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;}

static inline void no_data (Point point, scalar s) {int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;
  {foreach_child()
    val(s,0,0,0) = 1e30;end_foreach_child()}
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
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;







    return (27.*coarse(s,0,0,0) +
     9.*(coarse(s,child.x,0,0) + coarse(s,0,child.y,0) +
  coarse(s,0,0,child.z)) +
     3.*(coarse(s,child.x,child.y,0) + coarse(s,child.x,0,child.z) +
  coarse(s,0,child.y,child.z)) +
     coarse(s,child.x,child.y,child.z))/64.;

}

static inline void refine_bilinear (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;
  {foreach_child()
    val(s,0,0,0) = bilinear (point, s);end_foreach_child()}
}

static inline double quadratic (double a, double b, double c)
{
  return (30.*a + 5.*b - 3.*c)/32.;
}

static inline double biquadratic (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;
#line 169 "/home/lisergey/basilisk/src/grid/multigrid-common.h"
  if (!(false)) qassert ("/home/lisergey/basilisk/src/grid/multigrid-common.h", 169, "false");
  return 0.;

}

static inline double biquadratic_vertex (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;






  if (!(false)) qassert ("/home/lisergey/basilisk/src/grid/multigrid-common.h", 182, "false");
  return 0.;

}

static inline void refine_biquadratic (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;
  {foreach_child()
    val(s,0,0,0) = biquadratic (point, s);end_foreach_child()}
}

static inline void refine_linear (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;

#line 194
if(!is_constant(cm)){{
  coord g;
  if (_attribute[s.i].gradient)
    {
      g.x = _attribute[s.i].gradient (val(s,-1,0,0), val(s,0,0,0), val(s,1,0,0));
      
#line 198
g.y = _attribute[s.i].gradient (val(s,0,-1,0), val(s,0,0,0), val(s,0,1,0));
      
#line 198
g.z = _attribute[s.i].gradient (val(s,0,0,-1), val(s,0,0,0), val(s,0,0,1));}
  else
    {
      g.x = (val(s,1,0,0) - val(s,-1,0,0))/2.;
      
#line 201
g.y = (val(s,0,1,0) - val(s,0,-1,0))/2.;
      
#line 201
g.z = (val(s,0,0,1) - val(s,0,0,-1))/2.;}

  double sc = val(s,0,0,0), cmc = 4.*val(cm,0,0,0), sum = val(cm,0,0,0)*(1 << 3);
  {foreach_child() {
    val(s,0,0,0) = sc;
    
      val(s,0,0,0) += child.x*g.x*val(cm,-child.x,0,0)/cmc;
      
#line 207
val(s,0,0,0) += child.y*g.y*val(cm,0,-child.y,0)/cmc;
      
#line 207
val(s,0,0,0) += child.z*g.z*val(cm,0,0,-child.z)/cmc;
    sum -= val(cm,0,0,0);
  }end_foreach_child()}
  if (!(fabs(sum) < 1e-10)) qassert ("/home/lisergey/basilisk/src/grid/multigrid-common.h", 210, "fabs(sum) < 1e-10");
}}else {double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);

#line 194
{
  coord g;
  if (_attribute[s.i].gradient)
    {
      g.x = _attribute[s.i].gradient (val(s,-1,0,0), val(s,0,0,0), val(s,1,0,0));
      
#line 198
g.y = _attribute[s.i].gradient (val(s,0,-1,0), val(s,0,0,0), val(s,0,1,0));
      
#line 198
g.z = _attribute[s.i].gradient (val(s,0,0,-1), val(s,0,0,0), val(s,0,0,1));}
  else
    {
      g.x = (val(s,1,0,0) - val(s,-1,0,0))/2.;
      
#line 201
g.y = (val(s,0,1,0) - val(s,0,-1,0))/2.;
      
#line 201
g.z = (val(s,0,0,1) - val(s,0,0,-1))/2.;}

  double sc = val(s,0,0,0), cmc = 4.*_const_cm, sum = _const_cm*(1 << 3);
  {foreach_child() {
    val(s,0,0,0) = sc;
    
      val(s,0,0,0) += child.x*g.x*_const_cm/cmc;
      
#line 207
val(s,0,0,0) += child.y*g.y*_const_cm/cmc;
      
#line 207
val(s,0,0,0) += child.z*g.z*_const_cm/cmc;
    sum -= _const_cm;
  }end_foreach_child()}
  if (!(fabs(sum) < 1e-10)) qassert ("/home/lisergey/basilisk/src/grid/multigrid-common.h", 210, "fabs(sum) < 1e-10");
}}

#line 211
}

static inline void refine_reset (Point point, scalar v)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;
  {foreach_child()
    val(v,0,0,0) = 0.;end_foreach_child()}
}

static inline void refine_injection (Point point, scalar v)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;
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

static void multigrid_setup_vector (vector v)
{
   {
    _attribute[v.x.i].prolongation = refine_bilinear;
    _attribute[v.x.i].restriction = restriction_average;
  } 
#line 243
{
    _attribute[v.y.i].prolongation = refine_bilinear;
    _attribute[v.y.i].restriction = restriction_average;
  } 
#line 243
{
    _attribute[v.z.i].prolongation = refine_bilinear;
    _attribute[v.z.i].restriction = restriction_average;
  }
}

static vector multigrid_init_vector (vector v, const char * name)
{
  v = cartesian_init_vector (v, name);
  multigrid_setup_vector (v);
  return v;
}

static vector multigrid_init_face_vector (vector v, const char * name)
{
  v = cartesian_init_face_vector (v, name);
  
    _attribute[v.y.i].restriction = no_restriction;
    
#line 260
_attribute[v.z.i].restriction = no_restriction;
    
#line 260
_attribute[v.x.i].restriction = no_restriction;
  _attribute[v.x.i].restriction = restriction_face;
  return v;
}

static tensor multigrid_init_tensor (tensor t, const char * name)
{
  t = cartesian_init_tensor (t, name);
  
    multigrid_setup_vector (t.x);
    
#line 269
multigrid_setup_vector (t.y);
    
#line 269
multigrid_setup_vector (t.z);
  return t;
}

void multigrid_debug (Point point)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;
  cartesian_debug (point);

  FILE * plot = fopen ("plot", "a");
  if (point.level > 0) {
    char name[80] = "coarse";
    if (pid() > 0)
      sprintf (name, "coarse-%d", pid());
    FILE * fp = fopen (name, "w");
#line 307 "/home/lisergey/basilisk/src/grid/multigrid-common.h"
      double xc = x - child.x*Delta/2., yc = y - child.y*Delta/2.;
      double zc = z - child.z*Delta/2.;
      for (int k = 0; k <= 1; k++)
 for (int l = 0; l <= 1; l++)
   for (int m = 0; m <= 1; m++) {
     {scalar*_i=(scalar*)( all);if(_i)for(scalar v=*_i;(&v)->i>=0;v=*++_i){
       fprintf (fp, "%g %g %g %g ",
         xc + k*child.x*Delta*2. + _attribute[v.i].d.x*Delta,
         yc + l*child.y*Delta*2. + _attribute[v.i].d.y*Delta,
         zc + m*child.z*Delta*2. + _attribute[v.i].d.z*Delta,
         coarse(v,k*child.x,l*child.y,m*child.z));}}
     fputc ('\n', fp);
   }
      fprintf (ferr, ", '%s' u 1+4*v:2+4*v:3+4*v:4+4*v w labels tc lt 3 t ''",
        name);
      fprintf (plot, ", '%s' u 1+4*v:2+4*v:3+4*v:4+4*v w labels tc lt 3 t ''",
        name);

    fclose (fp);
  }

  if (is_coarse()) {
    char name[80] = "fine";
    if (pid() > 0)
      sprintf (name, "fine-%d", pid());
    FILE * fp = fopen (name, "w");
#line 364 "/home/lisergey/basilisk/src/grid/multigrid-common.h"
      double xf = x - Delta/4., yf = y - Delta/4., zf = z - Delta/4.;
      for (int k = -2; k <= 3; k++)
 for (int l = -2; l <= 3; l++)
   for (int m = -2; m <= 3; m++) {
     {scalar*_i=(scalar*)( all);if(_i)for(scalar v=*_i;(&v)->i>=0;v=*++_i){ {
       fprintf (fp, "%g %g %g ",
         xf + k*Delta/2. + _attribute[v.i].d.x*Delta/4.,
         yf + l*Delta/2. + _attribute[v.i].d.y*Delta/4.,
         zf + m*Delta/2. + _attribute[v.i].d.z*Delta/4.);
       if (allocated_child(k,l,m))
  fprintf (fp, "%g ", fine(v,k,l,m));
       else
  fputs ("n/a ", fp);
     }}}
     fputc ('\n', fp);
   }
      fprintf (ferr, ", '%s' u 1+4*v:2+4*v:3+4*v:4+4*v w labels tc lt 2 t ''",
        name);
      fprintf (plot, ", '%s' u 1+4*v:2+4*v:3+4*v:4+4*v w labels tc lt 2 t ''",
        name);

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
     
#line 404
list2 = list_add (list2, _attribute[s.i].v.y);
     
#line 404
list2 = list_add (list2, _attribute[s.i].v.z);}
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
    pfree (listdef,__func__,__FILE__,__LINE__);
    pfree (listc,__func__,__FILE__,__LINE__);
    pfree (list2,__func__,__FILE__,__LINE__);
  }
}

void multigrid_methods()
{
  cartesian_methods();
  init_scalar = multigrid_init_scalar;
  init_vertex_scalar = multigrid_init_vertex_scalar;
  init_vector = multigrid_init_vector;
  init_face_vector = multigrid_init_face_vector;
  init_tensor = multigrid_init_tensor;
  restriction = multigrid_restriction;
  debug = multigrid_debug;
}







void subtree_size (scalar size, bool leaves)
{




  foreach_stencil(1,{(NonLocal[]){{"size","scalar",(void *)&size,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n    \n// #line 454 \"/home/lisergey/basilisk/src/grid/multigrid-common.h\"\nval_out_(size,0,0,0) = 1;"})
    {_stencil_val_a(size,0,0,0);  }end_foreach_stencil();




  {
#line 453
foreach()
    val(size,0,0,0) = 1;end_foreach();}





  { Boundary ** _i = boundaries, * _b; while (_i && (_b = *_i++)) if (_b->restriction) _b->restriction (_b,((scalar[]) {size,{-1}}), depth()); };
  for (int l = depth() - 1; l >= 0; l--) {
    {foreach_coarse_level(l) {
      double sum = !leaves;
      {foreach_child()
 sum += val(size,0,0,0);end_foreach_child()}
      val(size,0,0,0) = sum;
    }end_foreach_coarse_level();}
    { Boundary ** _i = boundaries, * _b; while (_i && (_b = *_i++)) if (_b->restriction) _b->restriction (_b,((scalar[]) {size,{-1}}), l); };
  }
}
#line 5 "/home/lisergey/basilisk/src/grid/tree-common.h"




#line 21 "/home/lisergey/basilisk/src/grid/tree-common.h"
int refine_cell (Point point, scalar * list, int flag, Cache * refined)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;
  int nr = 0;


  if (level > 0)
    for (int k = 0; k != 2*child.x; k += child.x)

      for (int l = 0; l != 2*child.y; l += child.y)


 for (int m = 0; m != 2*child.z; m += child.z)

   if (aparent(k,l,m).pid >= 0 && is_leaf(aparent(k,l,m))) {
     Point p = point;


     p.level = point.level - 1;
     p.i = (point.i + 2)/2 + k;
     do { if (p.i < 2) p.i += 1 << p.level; else if (p.i >= 2 + (1 << p.level)) p.i -= 1 << p.level; } while(0);

       p.j = (point.j + 2)/2 + l;
       do { if (p.j < 2) p.j += 1 << p.level; else if (p.j >= 2 + (1 << p.level)) p.j -= 1 << p.level; } while(0);


       p.k = (point.k + 2)/2 + m;
       do { if (p.k < 2) p.k += 1 << p.level; else if (p.k >= 2 + (1 << p.level)) p.k -= 1 << p.level; } while(0);

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

#if 1
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
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;



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

#if 1
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
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;


  {foreach_child()
    if (cell.neighbors)
      {foreach_neighbor(1)
 if ((!is_leaf (cell) && cell.neighbors && cell.pid >= 0))
   coarsen_cell_recursive (point, list);end_foreach_neighbor()}end_foreach_child()}

  if (!(coarsen_cell (point, list))) qassert ("/home/lisergey/basilisk/src/grid/tree-common.h", 148, "coarsen_cell (point, list)");
}

void mpi_boundary_refine (scalar *);
void mpi_boundary_coarsen (int, int);
void mpi_boundary_update (scalar *);

static
scalar * list_add_depend (scalar * list, scalar s)
{
  if (is_constant(s) || _attribute[s.i].restriction == no_restriction)
    return list;
  {scalar*_i=(scalar*)( list);if(_i)for(scalar t=*_i;(&t)->i>=0;t=*++_i){
    if (t.i == s.i)
      return list;}}
  {scalar*_i=(scalar*)( _attribute[s.i].depends);if(_i)for(scalar d=*_i;(&d)->i>=0;d=*++_i){
    list = list_add_depend (list, d);}}
  return list_append (list, s);
}

typedef struct {
  int nc, nf;
} astats;

     
astats adapt_wavelet (scalar * slist,
        double * max,
        int maxlevel,
        int minlevel,
        scalar * list)
{tracing("adapt_wavelet","/home/lisergey/basilisk/src/grid/tree-common.h",173);
  scalar * ilist = list;

  if (is_constant(cm)) {
    if (list == NULL || list == all)
      list = list_copy (all);
    boundary_internal ((scalar *)list, "/home/lisergey/basilisk/src/grid/tree-common.h", 184);
    restriction (slist);
  }
  else {
    if (list == NULL || list == all) {
      list = list_copy (((scalar[]){cm, fm.x, fm.y, fm.z,{-1}}));
      {scalar*_i=(scalar*)( all);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
 list = list_add (list, s);}}
    }
    boundary_internal ((scalar *)list, "/home/lisergey/basilisk/src/grid/tree-common.h", 193);
    scalar * listr = list_concat (slist,((scalar[]) {cm,{-1}}));
    restriction (listr);
    pfree (listr,__func__,__FILE__,__LINE__);
  }

  astats st = {0, 0};
  scalar * listc = NULL;
  {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    listc = list_add_depend (listc, s);}}


  if (minlevel < 1)
    minlevel = 1;
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
   {scalar*_i=(scalar*)( slist);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {
     double emax = max[i++], sc[1 << 3];
     int c = 0;
     {foreach_child()
       sc[c++] = val(s,0,0,0);end_foreach_child()}
     _attribute[s.i].prolongation (point, s);
     c = 0;
     {foreach_child() {
       double e = fabs(sc[c] - val(s,0,0,0));
       if (e > emax && level < maxlevel) {
  cell.flags &= ~too_fine;
  cell.flags |= too_coarse;
       }
       else if ((e <= emax/1.5 || level > maxlevel) &&
         !(cell.flags & (too_coarse|just_fine))) {
  if (level >= minlevel)
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
       if (level >= maxlevel)
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
  pfree (listc,__func__,__FILE__,__LINE__);

  mpi_all_reduce (st.nf, MPI_INT, MPI_SUM);
  mpi_all_reduce (st.nc, MPI_INT, MPI_SUM);
  if (st.nc || st.nf)
    mpi_boundary_update (list);

  if (list != ilist)
    pfree (list,__func__,__FILE__,__LINE__);

  {end_tracing("adapt_wavelet","/home/lisergey/basilisk/src/grid/tree-common.h",316);return st;}
end_tracing("adapt_wavelet","/home/lisergey/basilisk/src/grid/tree-common.h",317);}
#line 339 "/home/lisergey/basilisk/src/grid/tree-common.h"
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
#line 384 "/home/lisergey/basilisk/src/grid/tree-common.h"
     
static void halo_face (vectorl vl)
{tracing("halo_face","/home/lisergey/basilisk/src/grid/tree-common.h",385);
  
    {scalar*_i=(scalar*)( vl.x);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      _attribute[s.i].dirty = 2;}}
    
#line 388
{scalar*_i=(scalar*)( vl.y);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      _attribute[s.i].dirty = 2;}}
    
#line 388
{scalar*_i=(scalar*)( vl.z);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      _attribute[s.i].dirty = 2;}}

  for (int l = depth() - 1; l >= 0; l--)
    {foreach_halo (prolongation, l)
      {
        if (vl.x) {
#line 410 "/home/lisergey/basilisk/src/grid/tree-common.h"
   if ((!is_leaf (neighbor(-1,0,0)) && neighbor(-1,0,0).neighbors && neighbor(-1,0,0).pid >= 0))
     {scalar*_i=(scalar*)( vl.x);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
       val(s,0,0,0) = (fine(s,0,0,0) + fine(s,0,1,0) +
       fine(s,0,0,1) + fine(s,0,1,1))/4.;}}
   if ((!is_leaf (neighbor(1,0,0)) && neighbor(1,0,0).neighbors && neighbor(1,0,0).pid >= 0))
     {scalar*_i=(scalar*)( vl.x);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
       val(s,1,0,0) = (fine(s,2,0,0) + fine(s,2,1,0) +
        fine(s,2,0,1) + fine(s,2,1,1))/4.;}}

 }
        
#line 394
if (vl.y) {
#line 410 "/home/lisergey/basilisk/src/grid/tree-common.h"
   if ((!is_leaf (neighbor(0,-1,0)) && neighbor(0,-1,0).neighbors && neighbor(0,-1,0).pid >= 0))
     {scalar*_i=(scalar*)( vl.y);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
       val(s,0,0,0) = (fine(s,0,0,0) + fine(s,0,0,1) +
       fine(s,1,0,0) + fine(s,1,0,1))/4.;}}
   if ((!is_leaf (neighbor(0,1,0)) && neighbor(0,1,0).neighbors && neighbor(0,1,0).pid >= 0))
     {scalar*_i=(scalar*)( vl.y);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
       val(s,0,1,0) = (fine(s,0,2,0) + fine(s,0,2,1) +
        fine(s,1,2,0) + fine(s,1,2,1))/4.;}}

 }
        
#line 394
if (vl.z) {
#line 410 "/home/lisergey/basilisk/src/grid/tree-common.h"
   if ((!is_leaf (neighbor(0,0,-1)) && neighbor(0,0,-1).neighbors && neighbor(0,0,-1).pid >= 0))
     {scalar*_i=(scalar*)( vl.z);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
       val(s,0,0,0) = (fine(s,0,0,0) + fine(s,1,0,0) +
       fine(s,0,1,0) + fine(s,1,1,0))/4.;}}
   if ((!is_leaf (neighbor(0,0,1)) && neighbor(0,0,1).neighbors && neighbor(0,0,1).pid >= 0))
     {scalar*_i=(scalar*)( vl.z);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
       val(s,0,0,1) = (fine(s,0,0,2) + fine(s,1,0,2) +
        fine(s,0,1,2) + fine(s,1,1,2))/4.;}}

 }}end_foreach_halo();}
end_tracing("halo_face","/home/lisergey/basilisk/src/grid/tree-common.h",420);}



static scalar tree_init_scalar (scalar s, const char * name)
{
  s = multigrid_init_scalar (s, name);
  _attribute[s.i].refine = _attribute[s.i].prolongation;
  return s;
}

static void prolongation_vertex (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;



  fine(s,1,1,1) = (val(s,0,0,0) + val(s,1,0,0) + val(s,0,1,0) + val(s,1,1,0) +
     val(s,0,0,1) + val(s,1,0,1) + val(s,0,1,1) + val(s,1,1,1))/8.;


  for (int i = 0; i <= 1; i++) {
    for (int j = 0; j <= 1; j++)

      for (int k = 0; k <= 1; k++)
 if (allocated_child(2*i,2*j,2*k))
   fine(s,2*i,2*j,2*k) = val(s,i,j,k);





    
      if (neighbor(i,0,0).neighbors) {



 fine(s,2*i,1,1) = (val(s,i,0,0) + val(s,i,1,0) + val(s,i,0,1) + val(s,i,1,1))/4.;
 fine(s,2*i,1,0) = (val(s,i,0,0) + val(s,i,1,0))/2.;
 fine(s,2*i,0,1) = (val(s,i,0,0) + val(s,i,0,1))/2.;
 if (allocated_child(2*i,1,2))
   fine(s,2*i,1,2) = (val(s,i,0,1) + val(s,i,1,1))/2.;
 if (allocated_child(2*i,2,1))
   fine(s,2*i,2,1) = (val(s,i,1,0) + val(s,i,1,1))/2.;

      }
      
#line 452
if (neighbor(0,i,0).neighbors) {



 fine(s,1,2*i,1) = (val(s,0,i,0) + val(s,0,i,1) + val(s,1,i,0) + val(s,1,i,1))/4.;
 fine(s,0,2*i,1) = (val(s,0,i,0) + val(s,0,i,1))/2.;
 fine(s,1,2*i,0) = (val(s,0,i,0) + val(s,1,i,0))/2.;
 if (allocated_child(2,2*i,1))
   fine(s,2,2*i,1) = (val(s,1,i,0) + val(s,1,i,1))/2.;
 if (allocated_child(1,2*i,2))
   fine(s,1,2*i,2) = (val(s,0,i,1) + val(s,1,i,1))/2.;

      }
      
#line 452
if (neighbor(0,0,i).neighbors) {



 fine(s,1,1,2*i) = (val(s,0,0,i) + val(s,1,0,i) + val(s,0,1,i) + val(s,1,1,i))/4.;
 fine(s,1,0,2*i) = (val(s,0,0,i) + val(s,1,0,i))/2.;
 fine(s,0,1,2*i) = (val(s,0,0,i) + val(s,0,1,i))/2.;
 if (allocated_child(1,2,2*i))
   fine(s,1,2,2*i) = (val(s,0,1,i) + val(s,1,1,i))/2.;
 if (allocated_child(2,1,2*i))
   fine(s,2,1,2*i) = (val(s,1,0,i) + val(s,1,1,i))/2.;

      }
  }
}

static scalar tree_init_vertex_scalar (scalar s, const char * name)
{
  s = multigrid_init_vertex_scalar (s, name);
  _attribute[s.i].refine = _attribute[s.i].prolongation = prolongation_vertex;
  return s;
}

static void tree_setup_vector (vector v)
{
  
    _attribute[v.x.i].refine = _attribute[v.x.i].prolongation;
    
#line 478
_attribute[v.y.i].refine = _attribute[v.y.i].prolongation;
    
#line 478
_attribute[v.z.i].refine = _attribute[v.z.i].prolongation;
}

static vector tree_init_vector (vector v, const char * name)
{
  v = multigrid_init_vector (v, name);
  tree_setup_vector (v);
  return v;
}


static void refine_face_x (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;
  vector v = _attribute[s.i].v;
#line 511 "/home/lisergey/basilisk/src/grid/tree-common.h"
  if (!(!is_leaf (neighbor(-1,0,0)) && neighbor(-1,0,0).neighbors && neighbor(-1,0,0).pid >= 0) &&
      (is_local(cell) || is_local(neighbor(-1,0,0)))) {
    double g1 = (val(v.x,0,+1,0) - val(v.x,0,-1,0))/8.;
    double g2 = (val(v.x,0,0,+1) - val(v.x,0,0,-1))/8.;
    for (int j = 0; j <= 1; j++)
      for (int k = 0; k <= 1; k++)
 fine(v.x,0,j,k) = val(v.x,0,0,0) + (2*j - 1)*g1 + (2*k - 1)*g2;
  }
  if (!(!is_leaf (neighbor(1,0,0)) && neighbor(1,0,0).neighbors && neighbor(1,0,0).pid >= 0) && neighbor(1,0,0).neighbors &&
      (is_local(cell) || is_local(neighbor(1,0,0)))) {
    double g1 = (val(v.x,1,+1,0) - val(v.x,1,-1,0))/8.;
    double g2 = (val(v.x,1,0,+1) - val(v.x,1,0,-1))/8.;
    for (int j = 0; j <= 1; j++)
      for (int k = 0; k <= 1; k++)
 fine(v.x,2,j,k) = val(v.x,1,0,0) + (2*j - 1)*g1 + (2*k - 1)*g2;
  }
  if (is_local(cell)) {
    double g1 = (val(v.x,0,+1,0) - val(v.x,0,-1,0) + val(v.x,1,+1,0) - val(v.x,1,-1,0))/16.;
    double g2 = (val(v.x,0,0,+1) - val(v.x,0,0,-1) + val(v.x,1,0,+1) - val(v.x,1,0,-1))/16.;
    for (int j = 0; j <= 1; j++)
      for (int k = 0; k <= 1; k++)
 fine(v.x,1,j,k) = (val(v.x,0,0,0) + val(v.x,1,0,0))/2. + (2*j - 1)*g1 + (2*k - 1)*g2;
  }

}

#line 489
static void refine_face_y (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;
  vector v = _attribute[s.i].v;
#line 511 "/home/lisergey/basilisk/src/grid/tree-common.h"
  if (!(!is_leaf (neighbor(0,-1,0)) && neighbor(0,-1,0).neighbors && neighbor(0,-1,0).pid >= 0) &&
      (is_local(cell) || is_local(neighbor(0,-1,0)))) {
    double g1 = (val(v.y,0,0,+1) - val(v.y,0,0,-1))/8.;
    double g2 = (val(v.y,+1,0,0) - val(v.y,-1,0,0))/8.;
    for (int j = 0; j <= 1; j++)
      for (int k = 0; k <= 1; k++)
 fine(v.y,k,0,j) = val(v.y,0,0,0) + (2*j - 1)*g1 + (2*k - 1)*g2;
  }
  if (!(!is_leaf (neighbor(0,1,0)) && neighbor(0,1,0).neighbors && neighbor(0,1,0).pid >= 0) && neighbor(0,1,0).neighbors &&
      (is_local(cell) || is_local(neighbor(0,1,0)))) {
    double g1 = (val(v.y,0,1,+1) - val(v.y,0,1,-1))/8.;
    double g2 = (val(v.y,+1,1,0) - val(v.y,-1,1,0))/8.;
    for (int j = 0; j <= 1; j++)
      for (int k = 0; k <= 1; k++)
 fine(v.y,k,2,j) = val(v.y,0,1,0) + (2*j - 1)*g1 + (2*k - 1)*g2;
  }
  if (is_local(cell)) {
    double g1 = (val(v.y,0,0,+1) - val(v.y,0,0,-1) + val(v.y,0,1,+1) - val(v.y,0,1,-1))/16.;
    double g2 = (val(v.y,+1,0,0) - val(v.y,-1,0,0) + val(v.y,+1,1,0) - val(v.y,-1,1,0))/16.;
    for (int j = 0; j <= 1; j++)
      for (int k = 0; k <= 1; k++)
 fine(v.y,k,1,j) = (val(v.y,0,0,0) + val(v.y,0,1,0))/2. + (2*j - 1)*g1 + (2*k - 1)*g2;
  }

}

#line 489
static void refine_face_z (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;
  vector v = _attribute[s.i].v;
#line 511 "/home/lisergey/basilisk/src/grid/tree-common.h"
  if (!(!is_leaf (neighbor(0,0,-1)) && neighbor(0,0,-1).neighbors && neighbor(0,0,-1).pid >= 0) &&
      (is_local(cell) || is_local(neighbor(0,0,-1)))) {
    double g1 = (val(v.z,+1,0,0) - val(v.z,-1,0,0))/8.;
    double g2 = (val(v.z,0,+1,0) - val(v.z,0,-1,0))/8.;
    for (int j = 0; j <= 1; j++)
      for (int k = 0; k <= 1; k++)
 fine(v.z,j,k,0) = val(v.z,0,0,0) + (2*j - 1)*g1 + (2*k - 1)*g2;
  }
  if (!(!is_leaf (neighbor(0,0,1)) && neighbor(0,0,1).neighbors && neighbor(0,0,1).pid >= 0) && neighbor(0,0,1).neighbors &&
      (is_local(cell) || is_local(neighbor(0,0,1)))) {
    double g1 = (val(v.z,+1,0,1) - val(v.z,-1,0,1))/8.;
    double g2 = (val(v.z,0,+1,1) - val(v.z,0,-1,1))/8.;
    for (int j = 0; j <= 1; j++)
      for (int k = 0; k <= 1; k++)
 fine(v.z,j,k,2) = val(v.z,0,0,1) + (2*j - 1)*g1 + (2*k - 1)*g2;
  }
  if (is_local(cell)) {
    double g1 = (val(v.z,+1,0,0) - val(v.z,-1,0,0) + val(v.z,+1,0,1) - val(v.z,-1,0,1))/16.;
    double g2 = (val(v.z,0,+1,0) - val(v.z,0,-1,0) + val(v.z,0,+1,1) - val(v.z,0,-1,1))/16.;
    for (int j = 0; j <= 1; j++)
      for (int k = 0; k <= 1; k++)
 fine(v.z,j,k,1) = (val(v.z,0,0,0) + val(v.z,0,0,1))/2. + (2*j - 1)*g1 + (2*k - 1)*g2;
  }

}

void refine_face (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;
  vector v = _attribute[s.i].v;
  
    _attribute[v.x.i].prolongation (point, v.x);
    
#line 541
_attribute[v.y.i].prolongation (point, v.y);
    
#line 541
_attribute[v.z.i].prolongation (point, v.z);
}

void refine_face_solenoidal (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;
  refine_face (point, s);

  if (is_local(cell)) {

    vector v = _attribute[s.i].v;
    double d[1 << 3], p[1 << 3];
    int i = 0;
    {foreach_child() {
      d[i] = 0.;
      
 d[i] += val(v.x,1,0,0) - val(v.x,0,0,0);
 
#line 556
d[i] += val(v.y,0,1,0) - val(v.y,0,0,0);
 
#line 556
d[i] += val(v.z,0,0,1) - val(v.z,0,0,0);
      i++;
    }end_foreach_child()}
#line 569 "/home/lisergey/basilisk/src/grid/tree-common.h"
    static double m[7][7] = {
      {7./12,5./24,3./8,5./24,3./8,1./4,1./3},
      {5./24,7./12,3./8,5./24,1./4,3./8,1./3},
      {3./8,3./8,3./4,1./4,3./8,3./8,1./2},
      {5./24,5./24,1./4,7./12,3./8,3./8,1./3},
      {3./8,1./4,3./8,3./8,3./4,3./8,1./2},
      {1./4,3./8,3./8,3./8,3./8,3./4,1./2},
      {1./3,1./3,1./2,1./3,1./2,1./2,5./6}
    };
    p[0] = 0.;
    for (int i = 0; i < 7; i++) {
      p[i + 1] = 0.;
      for (int j = 0; j < 7; j++)
 p[i + 1] += m[i][j]*d[j+1];
    }
    for (int k = 0; k <= 1; k++) {
      fine(v.x,1,0,k) += p[4+k] - p[0+k];
      fine(v.x,1,1,k) += p[6+k] - p[2+k];
      fine(v.y,0,1,k) += p[2+k] - p[0+k];
      fine(v.y,1,1,k) += p[6+k] - p[4+k];
    }
    fine(v.z,0,0,1) += p[1] - p[0];
    fine(v.z,0,1,1) += p[3] - p[2];
    fine(v.z,1,0,1) += p[5] - p[4];
    fine(v.z,1,1,1) += p[7] - p[6];

  }

}

static vector tree_init_face_vector (vector v, const char * name)
{
  v = cartesian_init_face_vector (v, name);
  
    _attribute[v.x.i].restriction = _attribute[v.x.i].refine = no_restriction;
    
#line 603
_attribute[v.y.i].restriction = _attribute[v.y.i].refine = no_restriction;
    
#line 603
_attribute[v.z.i].restriction = _attribute[v.z.i].refine = no_restriction;
  _attribute[v.x.i].restriction = restriction_face;
  _attribute[v.x.i].refine = refine_face;
  
    _attribute[v.x.i].prolongation = refine_face_x;
    
#line 607
_attribute[v.y.i].prolongation = refine_face_y;
    
#line 607
_attribute[v.z.i].prolongation = refine_face_z;
  return v;
}

static tensor tree_init_tensor (tensor t, const char * name)
{
  t = multigrid_init_tensor (t, name);
  
    tree_setup_vector (t.x);
    
#line 615
tree_setup_vector (t.y);
    
#line 615
tree_setup_vector (t.z);
  return t;
}

     
static void tree_boundary_level (scalar * list, int l)
{tracing("tree_boundary_level","/home/lisergey/basilisk/src/grid/tree-common.h",620);
  int depth = l < 0 ? depth() : l;

  if (tree_is_full()) {
    { Boundary ** _i = boundaries, * _b; while (_i && (_b = *_i++)) if (_b->level) _b->level (_b, list, depth); };
    {end_tracing("tree_boundary_level","/home/lisergey/basilisk/src/grid/tree-common.h",626);return;}
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
     
#line 640
list2 = list_add (list2, _attribute[s.i].v.y);
     
#line 640
list2 = list_add (list2, _attribute[s.i].v.z);}
 else {
   list2 = list_add (list2, s);
   if (_attribute[s.i].restriction == restriction_vertex)
     vlist = list_add (vlist, s);
 }
      }
    }}}

  if (vlist)
#line 674 "/home/lisergey/basilisk/src/grid/tree-common.h"
    {foreach_vertex () {
      if ((!is_leaf (cell) && cell.neighbors && cell.pid >= 0) || (!is_leaf (neighbor(-1,0,0)) && neighbor(-1,0,0).neighbors && neighbor(-1,0,0).pid >= 0) ||
   (!is_leaf (neighbor(0,-1,0)) && neighbor(0,-1,0).neighbors && neighbor(0,-1,0).pid >= 0) || (!is_leaf (neighbor(-1,-1,0)) && neighbor(-1,-1,0).neighbors && neighbor(-1,-1,0).pid >= 0) ||
   (!is_leaf (neighbor(0,0,-1)) && neighbor(0,0,-1).neighbors && neighbor(0,0,-1).pid >= 0) || (!is_leaf (neighbor(-1,0,-1)) && neighbor(-1,0,-1).neighbors && neighbor(-1,0,-1).pid >= 0) ||
   (!is_leaf (neighbor(0,-1,-1)) && neighbor(0,-1,-1).neighbors && neighbor(0,-1,-1).pid >= 0) || (!is_leaf (neighbor(-1,-1,-1)) && neighbor(-1,-1,-1).neighbors && neighbor(-1,-1,-1).pid >= 0)) {

 {scalar*_i=(scalar*)( vlist);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
   val(s,0,0,0) = is_vertex (child(0,0,0)) ? fine(s,0,0,0) : 1e30;}}
      }
      else
 { {
   if (child.y == 1 && child.z == 1 &&
       ((!is_leaf(cell) && !cell.neighbors && cell.pid >= 0) || (!is_leaf(neighbor(-1,0,0)) && !neighbor(-1,0,0).neighbors && neighbor(-1,0,0).pid >= 0))) {

     {scalar*_i=(scalar*)( vlist);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
       val(s,0,0,0) = is_vertex(neighbor(0,-1,-1)) && is_vertex(neighbor(0,1,-1))
  && is_vertex(neighbor(0,-1,1)) && is_vertex(neighbor(0,1,1)) ?
  (val(s,0,-1,-1) + val(s,0,1,-1) + val(s,0,-1,1) + val(s,0,1,1))/4. : 1e30;}}
   }
   else if (child.x == -1 && child.z == -1 && child.y == 1 &&
     ((!is_leaf(cell) && !cell.neighbors && cell.pid >= 0) || (!is_leaf(neighbor(-1,0,0)) && !neighbor(-1,0,0).neighbors && neighbor(-1,0,0).pid >= 0) ||
      (!is_leaf(neighbor(0,0,-1)) && !neighbor(0,0,-1).neighbors && neighbor(0,0,-1).pid >= 0) ||
      (!is_leaf(neighbor(-1,0,-1)) && !neighbor(-1,0,-1).neighbors && neighbor(-1,0,-1).pid >= 0))) {

     {scalar*_i=(scalar*)( vlist);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
       val(s,0,0,0) = is_vertex(neighbor(0,-1,0)) && is_vertex(neighbor(0,1,0)) ?
  (val(s,0,-1,0) + val(s,0,1,0))/2. : 1e30;}}
   }
 } 
#line 684
{
   if (child.z == 1 && child.x == 1 &&
       ((!is_leaf(cell) && !cell.neighbors && cell.pid >= 0) || (!is_leaf(neighbor(0,-1,0)) && !neighbor(0,-1,0).neighbors && neighbor(0,-1,0).pid >= 0))) {

     {scalar*_i=(scalar*)( vlist);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
       val(s,0,0,0) = is_vertex(neighbor(-1,0,-1)) && is_vertex(neighbor(-1,0,1))
  && is_vertex(neighbor(1,0,-1)) && is_vertex(neighbor(1,0,1)) ?
  (val(s,-1,0,-1) + val(s,-1,0,1) + val(s,1,0,-1) + val(s,1,0,1))/4. : 1e30;}}
   }
   else if (child.y == -1 && child.x == -1 && child.z == 1 &&
     ((!is_leaf(cell) && !cell.neighbors && cell.pid >= 0) || (!is_leaf(neighbor(0,-1,0)) && !neighbor(0,-1,0).neighbors && neighbor(0,-1,0).pid >= 0) ||
      (!is_leaf(neighbor(-1,0,0)) && !neighbor(-1,0,0).neighbors && neighbor(-1,0,0).pid >= 0) ||
      (!is_leaf(neighbor(-1,-1,0)) && !neighbor(-1,-1,0).neighbors && neighbor(-1,-1,0).pid >= 0))) {

     {scalar*_i=(scalar*)( vlist);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
       val(s,0,0,0) = is_vertex(neighbor(0,0,-1)) && is_vertex(neighbor(0,0,1)) ?
  (val(s,0,0,-1) + val(s,0,0,1))/2. : 1e30;}}
   }
 } 
#line 684
{
   if (child.x == 1 && child.y == 1 &&
       ((!is_leaf(cell) && !cell.neighbors && cell.pid >= 0) || (!is_leaf(neighbor(0,0,-1)) && !neighbor(0,0,-1).neighbors && neighbor(0,0,-1).pid >= 0))) {

     {scalar*_i=(scalar*)( vlist);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
       val(s,0,0,0) = is_vertex(neighbor(-1,-1,0)) && is_vertex(neighbor(1,-1,0))
  && is_vertex(neighbor(-1,1,0)) && is_vertex(neighbor(1,1,0)) ?
  (val(s,-1,-1,0) + val(s,1,-1,0) + val(s,-1,1,0) + val(s,1,1,0))/4. : 1e30;}}
   }
   else if (child.z == -1 && child.y == -1 && child.x == 1 &&
     ((!is_leaf(cell) && !cell.neighbors && cell.pid >= 0) || (!is_leaf(neighbor(0,0,-1)) && !neighbor(0,0,-1).neighbors && neighbor(0,0,-1).pid >= 0) ||
      (!is_leaf(neighbor(0,-1,0)) && !neighbor(0,-1,0).neighbors && neighbor(0,-1,0).pid >= 0) ||
      (!is_leaf(neighbor(0,-1,-1)) && !neighbor(0,-1,-1).neighbors && neighbor(0,-1,-1).pid >= 0))) {

     {scalar*_i=(scalar*)( vlist);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
       val(s,0,0,0) = is_vertex(neighbor(-1,0,0)) && is_vertex(neighbor(1,0,0)) ?
  (val(s,-1,0,0) + val(s,1,0,0))/2. : 1e30;}}
   }
 }}
    }end_foreach_vertex();}

  pfree (vlist,__func__,__FILE__,__LINE__);

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
    pfree (listdef,__func__,__FILE__,__LINE__);
    pfree (listc,__func__,__FILE__,__LINE__);
    pfree (list2,__func__,__FILE__,__LINE__);
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
     
#line 741
_attribute[v.y.i].prolongation (point, v.y);
     
#line 741
_attribute[v.z.i].prolongation (point, v.z);}}}
      }end_foreach_halo();}
      { Boundary ** _i = boundaries, * _b; while (_i && (_b = *_i++)) if (_b->level) _b->level (_b, list, i + 1); };
    }
    pfree (listr,__func__,__FILE__,__LINE__);
    pfree (listf,__func__,__FILE__,__LINE__);
  }
end_tracing("tree_boundary_level","/home/lisergey/basilisk/src/grid/tree-common.h",748);}

double treex (Point point) {int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;
  if (level == 0)
    return 0;




  if (!(false)) qassert ("/home/lisergey/basilisk/src/grid/tree-common.h", 757, "false");
  double i = 0;

  return treex(parent) + i/(1 << 2*(level - 1));
}

double treey (Point point) {int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;
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
{tracing("tree_check","/home/lisergey/basilisk/src/grid/tree-common.h",780);


  long nleaves = 0, nactive = 0;
  {foreach_cell_all() {
    if (is_leaf(cell)) {
      if (!(cell.pid >= 0)) qassert ("/home/lisergey/basilisk/src/grid/tree-common.h", 787, "cell.pid >= 0");
      nleaves++;
    }
    if (is_local(cell))
      if (!(is_active(cell) || (!is_leaf(cell) && !cell.neighbors && cell.pid >= 0))) qassert ("/home/lisergey/basilisk/src/grid/tree-common.h", 791, "is_active(cell) || is_prolongation(cell)");
    if (is_active(cell))
      nactive++;

    int neighbors = 0;
    {foreach_neighbor(1)
      if (allocated(0,0,0) && (!is_leaf (cell) && cell.neighbors && cell.pid >= 0))
 neighbors++;end_foreach_neighbor()}
    if (!(cell.neighbors == neighbors)) qassert ("/home/lisergey/basilisk/src/grid/tree-common.h", 799, "cell.neighbors == neighbors");


    if (!cell.neighbors)
      if (!(!allocated_child(0,0,0))) qassert ("/home/lisergey/basilisk/src/grid/tree-common.h", 803, "!allocated_child(0)");
  }end_foreach_cell_all();}


  long reachable = 0;
  {foreach_cell() {
    if (is_active(cell))
      reachable++;
    else
      continue;
  }end_foreach_cell();}
  if (!(nactive == reachable)) qassert ("/home/lisergey/basilisk/src/grid/tree-common.h", 814, "nactive == reachable");


  reachable = 0;
  {foreach_cell()
    if (is_leaf(cell)) {
      reachable++;
      continue;
    }end_foreach_cell();}
  if (!(nleaves == reachable)) qassert ("/home/lisergey/basilisk/src/grid/tree-common.h", 823, "nleaves == reachable");
end_tracing("tree_check","/home/lisergey/basilisk/src/grid/tree-common.h",824);}

     
static void tree_restriction (scalar * list) {tracing("tree_restriction","/home/lisergey/basilisk/src/grid/tree-common.h",827);
  boundary_internal ((scalar *)list, "/home/lisergey/basilisk/src/grid/tree-common.h", 828);
  if (tree_is_full())
    multigrid_restriction (list);
end_tracing("tree_restriction","/home/lisergey/basilisk/src/grid/tree-common.h",831);}

void tree_methods()
{
  multigrid_methods();
  init_scalar = tree_init_scalar;
  init_vertex_scalar = tree_init_vertex_scalar;
  init_vector = tree_init_vector;
  init_face_vector = tree_init_face_vector;
  init_tensor = tree_init_tensor;
  boundary_level = tree_boundary_level;
  boundary_face = halo_face;
  restriction = tree_restriction;
}
#line 1670 "/home/lisergey/basilisk/src/grid/tree.h"


void tree_periodic (int dir)
{
  int depth = grid ? depth() : -1;
  if (grid)
    free_grid();
  periodic (dir);
  if (depth >= 0)
    init_grid (1 << depth);
}


#if 1
#line 1 "grid/tree-mpi.h"
#line 1 "/home/lisergey/basilisk/src/grid/tree-mpi.h"

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
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;
  if (level > rcv->depth) {
    rcv->halo = (CacheLevel *) prealloc (rcv->halo, (level + 1)*sizeof(CacheLevel),__func__,__FILE__,__LINE__);
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
    pfree (rcv->buf,__func__,__FILE__,__LINE__);
    rcv->buf = NULL;
    prof_stop();
  }
}

static void rcv_destroy (Rcv * rcv)
{
  rcv_free_buf (rcv);
  for (int i = 0; i <= rcv->depth; i++)
    if (rcv->halo[i].n > 0)
      pfree (rcv->halo[i].p,__func__,__FILE__,__LINE__);
  pfree (rcv->halo,__func__,__FILE__,__LINE__);
}

static RcvPid * rcv_pid_new (const char * name)
{
  RcvPid * r = ((RcvPid *) pcalloc (1, sizeof(RcvPid),__func__,__FILE__,__LINE__));
  r->name = pstrdup (name,__func__,__FILE__,__LINE__);
  return r;
}

static Rcv * rcv_pid_pointer (RcvPid * p, int pid)
{
  if (!(pid >= 0 && pid < npe())) qassert ("/home/lisergey/basilisk/src/grid/tree-mpi.h", 88, "pid >= 0 && pid < npe()");

  int i;
  for (i = 0; i < p->npid; i++)
    if (pid == p->rcv[i].pid)
      break;

  if (i == p->npid) {
    p->rcv = (Rcv *) prealloc (p->rcv, (++p->npid)*sizeof(Rcv),__func__,__FILE__,__LINE__);
    Rcv * rcv = &p->rcv[p->npid-1];
    rcv->pid = pid;
    rcv->depth = rcv->maxdepth = 0;
    rcv->halo = ((CacheLevel *) pmalloc ((1)*sizeof(CacheLevel),__func__,__FILE__,__LINE__));
    rcv->buf = NULL;
    cache_level_init (&rcv->halo[0]);
  }
  return &p->rcv[i];
}

static void rcv_pid_append (RcvPid * p, int pid, Point point)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;
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
  pfree (p->rcv,__func__,__FILE__,__LINE__);
  pfree (p->name,__func__,__FILE__,__LINE__);
  pfree (p,__func__,__FILE__,__LINE__);
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
 if (*b != 1e30 && allocated(1,0,0))
   memcpy (&val(v.x,1,0,0), b, sizeof(double)*_attribute[v.x.i].block);
 b += _attribute[v.x.i].block;
      } 
#line 171
{
 memcpy (&val(v.y,0,0,0), b, sizeof(double)*_attribute[v.y.i].block);
 b += _attribute[v.y.i].block;
 if (*b != 1e30 && allocated(0,1,0))
   memcpy (&val(v.y,0,1,0), b, sizeof(double)*_attribute[v.y.i].block);
 b += _attribute[v.y.i].block;
      } 
#line 171
{
 memcpy (&val(v.z,0,0,0), b, sizeof(double)*_attribute[v.z.i].block);
 b += _attribute[v.z.i].block;
 if (*b != 1e30 && allocated(0,0,1))
   memcpy (&val(v.z,0,0,1), b, sizeof(double)*_attribute[v.z.i].block);
 b += _attribute[v.z.i].block;
      }}}}
    {scalar*_i=(scalar*)( listv);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {
      for (int i = 0; i <= 1; i++)
 for (int j = 0; j <= 1; j++)

   for (int k = 0; k <= 1; k++) {
     if (*b != 1e30 && allocated(i,j,k))
       memcpy (&val(s,i,j,k), b, sizeof(double)*_attribute[s.i].block);
     b += _attribute[s.i].block;
   }







    }}}
  }end_foreach_cache_level();}
  size_t size = b - (double *) rcv->buf;
  pfree (rcv->buf,__func__,__FILE__,__LINE__);
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
#line 234 "/home/lisergey/basilisk/src/grid/tree-mpi.h"
static void mpi_recv_check (void * buf, int count, MPI_Datatype datatype,
       int source, int tag,
       MPI_Comm comm, MPI_Status * status,
       const char * name)
{
#line 269 "/home/lisergey/basilisk/src/grid/tree-mpi.h"
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
{tracing("mpi_waitany","/home/lisergey/basilisk/src/grid/tree-mpi.h",291);
  { int _ret= MPI_Waitany (count, array_of_requests, indx, status);end_tracing("mpi_waitany","/home/lisergey/basilisk/src/grid/tree-mpi.h",294);return _ret;}
end_tracing("mpi_waitany","/home/lisergey/basilisk/src/grid/tree-mpi.h",295);}

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

  int len = list_lenb (list) + 2*3*vectors_lenb (listf) +
    (1 << 3)*list_lenb (listv);

  MPI_Request r[m->npid];
  Rcv * rrcv[m->npid];
  int nr = 0;
  for (int i = 0; i < m->npid; i++) {
    Rcv * rcv = &m->rcv[i];
    if (l <= rcv->depth && rcv->halo[l].n > 0) {
      if (!(!rcv->buf)) qassert ("/home/lisergey/basilisk/src/grid/tree-mpi.h", 328, "!rcv->buf");
      rcv->buf = pmalloc (sizeof (double)*rcv->halo[l].n*len,__func__,__FILE__,__LINE__);






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
      if (!(l <= rcv->depth && rcv->halo[l].n > 0)) qassert ("/home/lisergey/basilisk/src/grid/tree-mpi.h", 355, "l <= rcv->depth && rcv->halo[l].n > 0");
      if (!(rcv->buf)) qassert ("/home/lisergey/basilisk/src/grid/tree-mpi.h", 356, "rcv->buf");
      apply_bc (rcv, list, listv, listf, l, s);
      mpi_waitany (nr, r, &i, &s);
    }
  }

  prof_stop();
}

     
static void rcv_pid_wait (RcvPid * m)
{tracing("rcv_pid_wait","/home/lisergey/basilisk/src/grid/tree-mpi.h",366);

  for (int i = 0; i < m->npid; i++)
    rcv_free_buf (&m->rcv[i]);
end_tracing("rcv_pid_wait","/home/lisergey/basilisk/src/grid/tree-mpi.h",371);}

static void rcv_pid_send (RcvPid * m, scalar * list, scalar * listv,
     vector * listf, int l)
{
  if (m->npid == 0)
    return;

  prof_start ("rcv_pid_send");

  int len = list_lenb (list) + 2*3*vectors_lenb (listf) +
    (1 << 3)*list_lenb (listv);


  for (int i = 0; i < m->npid; i++) {
    Rcv * rcv = &m->rcv[i];
    if (l <= rcv->depth && rcv->halo[l].n > 0) {
      if (!(!rcv->buf)) qassert ("/home/lisergey/basilisk/src/grid/tree-mpi.h", 388, "!rcv->buf");
      rcv->buf = pmalloc (sizeof (double)*rcv->halo[l].n*len,__func__,__FILE__,__LINE__);
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
       *b = 1e30;
     b += _attribute[v.x.i].block;
   } 
#line 397
{
     memcpy (b, &val(v.y,0,0,0), sizeof(double)*_attribute[v.y.i].block);
     b += _attribute[v.y.i].block;
     if (allocated(0,1,0))
       memcpy (b, &val(v.y,0,1,0), sizeof(double)*_attribute[v.y.i].block);
     else
       *b = 1e30;
     b += _attribute[v.y.i].block;
   } 
#line 397
{
     memcpy (b, &val(v.z,0,0,0), sizeof(double)*_attribute[v.z.i].block);
     b += _attribute[v.z.i].block;
     if (allocated(0,0,1))
       memcpy (b, &val(v.z,0,0,1), sizeof(double)*_attribute[v.z.i].block);
     else
       *b = 1e30;
     b += _attribute[v.z.i].block;
   }}}}
 {scalar*_i=(scalar*)( listv);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {
   for (int i = 0; i <= 1; i++)
     for (int j = 0; j <= 1; j++)

       for (int k = 0; k <= 1; k++) {
  if (allocated(i,j,k))
    memcpy (b, &val(s,i,j,k), sizeof(double)*_attribute[s.i].block);
  else
    *b = 1e30;
  b += _attribute[s.i].block;
       }
#line 426 "/home/lisergey/basilisk/src/grid/tree-mpi.h"
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
  pfree (listr,__func__,__FILE__,__LINE__);
  pfree (listf,__func__,__FILE__,__LINE__);
  pfree (listv,__func__,__FILE__,__LINE__);
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
  pfree (m,__func__,__FILE__,__LINE__);
}

     
static void mpi_boundary_level (const Boundary * b, scalar * list, int l)
{tracing("mpi_boundary_level","/home/lisergey/basilisk/src/grid/tree-mpi.h",492);
  MpiBoundary * m = (MpiBoundary *) b;
  rcv_pid_sync (&m->mpi_level, list, l);
  rcv_pid_sync (&m->mpi_level_root, list, l);
end_tracing("mpi_boundary_level","/home/lisergey/basilisk/src/grid/tree-mpi.h",497);}

     
static void mpi_boundary_restriction (const Boundary * b, scalar * list, int l)
{tracing("mpi_boundary_restriction","/home/lisergey/basilisk/src/grid/tree-mpi.h",500);
  MpiBoundary * m = (MpiBoundary *) b;
  rcv_pid_sync (&m->restriction, list, l);
end_tracing("mpi_boundary_restriction","/home/lisergey/basilisk/src/grid/tree-mpi.h",504);}

void mpi_boundary_new()
{
  mpi_boundary = (Boundary *) ((MpiBoundary *) pcalloc (1, sizeof(MpiBoundary),__func__,__FILE__,__LINE__));
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
    fprintf (fp, "%s%g %g %g %d\n", prefix, x, y, z, level);}end_is_face_y()
#line 581
is_face_z(){
    fprintf (fp, "%s%g %g %g %d\n", prefix, x, y, z, level);}end_is_face_z()}end_foreach_face_generic();}
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
    if (!(cell.neighbors == n)) qassert ("/home/lisergey/basilisk/src/grid/tree-mpi.h", 599, "cell.neighbors == n");
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
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;
  if ((!is_leaf (cell) && cell.neighbors && cell.pid >= 0))
    {foreach_child()
      if (is_local(cell))
 return true;end_foreach_child()}
  return false;
}


static bool is_local_prolongation (Point point, Point p)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;



  struct { int x, y, z; } dp = {p.i - point.i, p.j - point.j, p.k - point.k};

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
#line 713
{
    if (dp.z == 0 && ((!is_leaf (neighbor(0,0,-1)) && neighbor(0,0,-1).neighbors && neighbor(0,0,-1).pid >= 0) || (!is_leaf (neighbor(0,0,1)) && neighbor(0,0,1).neighbors && neighbor(0,0,1).pid >= 0)))
      return true;
    if ((!is_leaf (neighbor(0,0,dp.z)) && neighbor(0,0,dp.z).neighbors && neighbor(0,0,dp.z).pid >= 0))
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
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;
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
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;
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
  int * row = ((int *) pmalloc ((npe())*sizeof(int),__func__,__FILE__,__LINE__));
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

      int * astatus = ((int *) pmalloc ((npe())*sizeof(int),__func__,__FILE__,__LINE__));
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
      pfree (astatus,__func__,__FILE__,__LINE__);

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
  pfree (row,__func__,__FILE__,__LINE__);
}

static bool has_local_child (Point point)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;
  {foreach_child()
    if (is_local(cell))
      return true;end_foreach_child()}
  return false;
}

     
void mpi_boundary_update_buffers()
{tracing("mpi_boundary_update_buffers","/home/lisergey/basilisk/src/grid/tree-mpi.h",858);
  if (npe() == 1)
    {end_tracing("mpi_boundary_update_buffers","/home/lisergey/basilisk/src/grid/tree-mpi.h",861);return;}

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
 pfree (pids.p,__func__,__FILE__,__LINE__);
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
     pfree (pids.p,__func__,__FILE__,__LINE__);
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
   pfree (pids.p,__func__,__FILE__,__LINE__);
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
#line 1015 "/home/lisergey/basilisk/src/grid/tree-mpi.h"
end_tracing("mpi_boundary_update_buffers","/home/lisergey/basilisk/src/grid/tree-mpi.h",1015);}

     
void mpi_boundary_refine (scalar * list)
{tracing("mpi_boundary_refine","/home/lisergey/basilisk/src/grid/tree-mpi.h",1018);
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


  pfree (((Tree *)grid)->refined.p,__func__,__FILE__,__LINE__);
  ((Tree *)grid)->refined = rerefined;

  prof_stop();



  mpi_all_reduce (rerefined.n, MPI_INT, MPI_SUM);
  if (rerefined.n)
    mpi_boundary_refine (list);
  {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    _attribute[s.i].dirty = true;}}
end_tracing("mpi_boundary_refine","/home/lisergey/basilisk/src/grid/tree-mpi.h",1086);}

static void check_depth()
{
#line 1121 "/home/lisergey/basilisk/src/grid/tree-mpi.h"
}

typedef struct {
  int refined, leaf;
} Remote;



     
void mpi_boundary_coarsen (int l, int too_fine)
{tracing("mpi_boundary_coarsen","/home/lisergey/basilisk/src/grid/tree-mpi.h",1130);
  if (npe() == 1)
    {end_tracing("mpi_boundary_coarsen","/home/lisergey/basilisk/src/grid/tree-mpi.h",1133);return;}

  check_depth();

  if (!(sizeof(Remote) == sizeof(double))) qassert ("/home/lisergey/basilisk/src/grid/tree-mpi.h", 1137, "sizeof(Remote) == sizeof(double)");

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
  mpi_boundary_level (mpi_boundary,((scalar[]) {remote,{-1}}), l);

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
    mpi_boundary_level (mpi_boundary,((scalar[]) {remote,{-1}}), l);
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
end_tracing("mpi_boundary_coarsen","/home/lisergey/basilisk/src/grid/tree-mpi.h",1196);}

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
{tracing("mpi_partitioning","/home/lisergey/basilisk/src/grid/tree-mpi.h",1253);
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
end_tracing("mpi_partitioning","/home/lisergey/basilisk/src/grid/tree-mpi.h",1293);}

void restore_mpi (FILE * fp, scalar * list1)
{
  long index = 0, nt = 0, start = ftell (fp);
  scalar  size=new_scalar("size"), * list = list_concat (((scalar[]){size,{-1}}), list1);;
  long offset = sizeof(double)*list_len(list);


  static const unsigned short set = 1 << user;
  scalar * listm = is_constant(cm) ? NULL : (scalar *)((vector[]){fm,{{-1},{-1},{-1}}});
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
  pfree (list,__func__,__FILE__,__LINE__);delete((scalar*)((scalar[]){size,{-1}}));
}
#line 1435 "/home/lisergey/basilisk/src/grid/tree-mpi.h"
     
double z_indexing (scalar index, bool leaves)
{tracing("z_indexing","/home/lisergey/basilisk/src/grid/tree-mpi.h",1436);



  scalar  size=new_scalar("size");
  subtree_size (size, leaves);






  double maxi = -1.;
  if (pid() == 0)
    {foreach_level(0)
      maxi = val(size,0,0,0) - 1.;end_foreach_level();}




  {foreach_level(0)
    val(index,0,0,0) = 0;end_foreach_level();}
  for (int l = 0; l < depth(); l++) {
    { Boundary ** _i = boundaries, * _b; while (_i && (_b = *_i++)) if (_b->restriction) _b->restriction (_b,((scalar[]) {index,{-1}}), l); };
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
  { Boundary ** _i = boundaries, * _b; while (_i && (_b = *_i++)) if (_b->restriction) _b->restriction (_b,((scalar[]) {index,{-1}}), depth()); };

  {delete((scalar*)((scalar[]){size,{-1}}));{end_tracing("z_indexing","/home/lisergey/basilisk/src/grid/tree-mpi.h",1493);return maxi;}}delete((scalar*)((scalar[]){size,{-1}}));
end_tracing("z_indexing","/home/lisergey/basilisk/src/grid/tree-mpi.h",1494);}
#line 1685 "/home/lisergey/basilisk/src/grid/tree.h"
#line 1 "grid/balance.h"
#line 1 "/home/lisergey/basilisk/src/grid/balance.h"


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
 if (!(is_leaf(cell))) qassert ("/home/lisergey/basilisk/src/grid/balance.h", 34, "is_leaf(cell)");
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
      if (!(c->neighbors)) qassert ("/home/lisergey/basilisk/src/grid/balance.h", 81, "c->neighbors");\
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
    a.p = pmalloc (a.len,__func__,__FILE__,__LINE__);
    if (fp)
      fprintf (fp, "receiving %ld from %d\n", a.len, from);
    mpi_recv_check (a.p, a.len, MPI_BYTE, from, (256),
      MPI_COMM_WORLD, MPI_STATUS_IGNORE, "receive_tree (p)");

    {foreach_tree (&a, sizeof(Cell) + datasize, NULL) {
      memcpy (((char *)&cell) + sizeof(Cell), ((char *)c) + sizeof(Cell),
       datasize);
      if (!(((NewPid *)&val(newpid,0,0,0))->pid > 0)) qassert ("/home/lisergey/basilisk/src/grid/balance.h", 160, "NEWPID()->pid > 0");
      if (fp)
 fprintf (fp, "%g %g %g %d %d %d %d %d %d recv\n",
   x, y, z, ((NewPid *)&val(newpid,0,0,0))->pid - 1, cell.pid,
   c->flags & leaf,
   cell.flags & leaf, from, ((NewPid *)&val(newpid,0,0,0))->leaf);
    }end_foreach_tree();}
    pfree (a.p,__func__,__FILE__,__LINE__);
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
{tracing("balance","/home/lisergey/basilisk/src/grid/balance.h",201);
  if (npe() == 1)
    {end_tracing("balance","/home/lisergey/basilisk/src/grid/balance.h",204);return false;}

  if (!(sizeof(NewPid) == sizeof(double))) qassert ("/home/lisergey/basilisk/src/grid/balance.h", 206, "sizeof(NewPid) == sizeof(double)");

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
    {end_tracing("balance","/home/lisergey/basilisk/src/grid/balance.h",244);return false;}

  scalar  newpid=new_scalar("newpid");
  double zn = z_indexing (newpid, mpi.leaves);
  if (pid() == 0)
    if (!(zn + 1 == nt)) qassert ("/home/lisergey/basilisk/src/grid/balance.h", 249, "zn + 1 == nt");

  FILE * fp = NULL;
#line 261 "/home/lisergey/basilisk/src/grid/balance.h"
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
    { Boundary ** _i = boundaries, * _b; while (_i && (_b = *_i++)) if (_b->level) _b->level (_b,((scalar[]) {newpid,{-1}}), l); };
#line 305 "/home/lisergey/basilisk/src/grid/balance.h"
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

  {delete((scalar*)((scalar[]){newpid,{-1}}));{end_tracing("balance","/home/lisergey/basilisk/src/grid/balance.h",392);return pid_changed;}}delete((scalar*)((scalar[]){newpid,{-1}}));
end_tracing("balance","/home/lisergey/basilisk/src/grid/balance.h",393);}

void mpi_boundary_update (scalar * list)
{
  mpi_boundary_update_buffers();
  {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    _attribute[s.i].dirty = true;}}
  grid->tn = 0;
  boundary_internal ((scalar *)list, "/home/lisergey/basilisk/src/grid/balance.h", 401);
  while (balance());
}
#line 1686 "/home/lisergey/basilisk/src/grid/tree.h"
#else
void mpi_boundary_refine (scalar * list){}
void mpi_boundary_coarsen (int a, int b){}
void mpi_boundary_update (scalar * list) {
  {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    _attribute[s.i].dirty = true;}}
  boundary_internal ((scalar *)list, "/home/lisergey/basilisk/src/grid/tree.h", 1692);
}
#endif
#line 4 "/home/lisergey/basilisk/src/grid/octree.h"

void octree_methods() {
  tree_methods();
}
#line 8 "cylinder.c"
#line 1 "fractions.h"
#line 1 "/home/lisergey/basilisk/src/fractions.h"
#line 12 "/home/lisergey/basilisk/src/fractions.h"
#line 1 "geometry.h"
#line 1 "/home/lisergey/basilisk/src/geometry.h"
#line 35 "/home/lisergey/basilisk/src/geometry.h"
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



double plane_alpha (double c, coord n)
{
  double alpha;
  coord n1;

  n1.x = fabs (n.x); n1.y = fabs (n.y); n1.z = fabs (n.z);

  double m1, m2, m3;
  m1 = min(n1.x, n1.y);
  m3 = max(n1.x, n1.y);
  m2 = n1.z;
  if (m2 < m1) {
    double tmp = m1;
    m1 = m2;
    m2 = tmp;
  }
  else if (m2 > m3) {
    double tmp = m3;
    m3 = m2;
    m2 = tmp;
  }
  double m12 = m1 + m2;
  double pr = max(6.*m1*m2*m3, 1e-50);
  double V1 = m1*m1*m1/pr;
  double V2 = V1 + (m2 - m1)/(2.*m3), V3;
  double mm;
  if (m3 < m12) {
    mm = m3;
    V3 = (m3*m3*(3.*m12 - m3) + m1*m1*(m1 - 3.*m3) + m2*m2*(m2 - 3.*m3))/pr;
  }
  else {
    mm = m12;
    V3 = mm/(2.*m3);
  }

  c = clamp (c, 0., 1.);
  double ch = min(c, 1. - c);
  if (ch < V1)
    alpha = pow (pr*ch, 1./3.);
  else if (ch < V2)
    alpha = (m1 + sqrt(m1*m1 + 8.*m2*m3*(ch - V1)))/2.;
  else if (ch < V3) {
    double p12 = sqrt (2.*m1*m2);
    double q = 3.*(m12 - 2.*m3*ch)/(4.*p12);
    double teta = acos(clamp(q,-1.,1.))/3.;
    double cs = cos(teta);
    alpha = p12*(sqrt(3.*(1. - cs*cs)) - cs) + m12;
  }
  else if (m12 <= m3)
    alpha = m3*ch + mm/2.;
  else {
    double p = m1*(m2 + m3) + m2*m3 - 1./4., p12 = sqrt(p);
    double q = 3.*m1*m2*m3*(1./2. - ch)/(2.*p*p12);
    double teta = acos(clamp(q,-1.,1.))/3.;
    double cs = cos(teta);
    alpha = p12*(sqrt(3.*(1. - cs*cs)) - cs) + 1./2.;
  }
  if (c > 1./2.) alpha = 1. - alpha;

  if (n.x < 0.)
    alpha += n.x;
  if (n.y < 0.)
    alpha += n.y;
  if (n.z < 0.)
    alpha += n.z;

  return alpha - (n.x + n.y + n.z)/2.;;
}
#line 163 "/home/lisergey/basilisk/src/geometry.h"
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



double plane_volume (coord n, double alpha)
{
  double al = alpha + (n.x + n.y + n.z)/2. +
    max(0., -n.x) + max(0., -n.y) + max(0., -n.z);
  if (al <= 0.)
    return 0.;
  double tmp = fabs(n.x) + fabs(n.y) + fabs(n.z);
  if (al >= tmp)
    return 1.;
  if (tmp < 1e-10)
    return 0.;
  double n1 = fabs(n.x)/tmp;
  double n2 = fabs(n.y)/tmp;
  double n3 = fabs(n.z)/tmp;
  al = max(0., min(1., al/tmp));
  double al0 = min(al, 1. - al);
  double b1 = min(n1, n2);
  double b3 = max(n1, n2);
  double b2 = n3;
  if (b2 < b1) {
    tmp = b1;
    b1 = b2;
    b2 = tmp;
  }
  else if (b2 > b3) {
    tmp = b3;
    b3 = b2;
    b2 = tmp;
  }
  double b12 = b1 + b2;
  double bm = min(b12, b3);
  double pr = max(6.*b1*b2*b3, 1e-50);
  if (al0 < b1)
    tmp = al0*al0*al0/pr;
  else if (al0 < b2)
    tmp = 0.5*al0*(al0 - b1)/(b2*b3) + b1*b1*b1/pr;
  else if (al0 < bm)
    tmp = (al0*al0*(3.*b12 - al0) + b1*b1*(b1 - 3.*al0) +
    b2*b2*(b2 - 3.*al0))/pr;
  else if (b12 < b3)
    tmp = (al0 - 0.5*bm)/b3;
  else
    tmp = (al0*al0*(3. - 2.*al0) + b1*b1*(b1 - 3.*al0) +
    b2*b2*(b2 - 3.*al0) + b3*b3*(b3 - 3.*al0))/pr;

  double volume = al <= 0.5 ? tmp : 1. - tmp;
  return clamp (volume, 0., 1.);
}
#line 267 "/home/lisergey/basilisk/src/geometry.h"
double rectangle_fraction (coord n, double alpha, coord a, coord b)
{
  coord n1;
   {
    alpha -= n.x*(b.x + a.x)/2.;
    n1.x = n.x*(b.x - a.x);
  } 
#line 270
{
    alpha -= n.y*(b.y + a.y)/2.;
    n1.y = n.y*(b.y - a.y);
  } 
#line 270
{
    alpha -= n.z*(b.z + a.z)/2.;
    n1.z = n.z*(b.z - a.z);
  }
  return plane_volume (n1, alpha);
}
#line 307 "/home/lisergey/basilisk/src/geometry.h"
static coord cube_edge[12][2] = {
  {{0.,0.,0.},{1.,0.,0.}},{{0.,0.,1.},{1.,0.,1.}},
  {{0.,1.,1.},{1.,1.,1.}},{{0.,1.,0.},{1.,1.,0.}},
  {{0.,0.,0.},{0.,1.,0.}},{{0.,0.,1.},{0.,1.,1.}},
  {{1.,0.,1.},{1.,1.,1.}},{{1.,0.,0.},{1.,1.,0.}},
  {{0.,0.,0.},{0.,0.,1.}},{{1.,0.,0.},{1.,0.,1.}},
  {{1.,1.,0.},{1.,1.,1.}},{{0.,1.,0.},{0.,1.,1.}}
};




static int cube_connect[12][2][4] = {
  {{9, 1, 8}, {4, 3, 7}},
  {{6, 2, 5}, {8, 0, 9}},
  {{10, 3, 11}, {5, 1, 6}},
  {{7, 0, 4}, {11, 2, 10}},
  {{3, 7, 0}, {8, 5, 11}},
  {{11, 4, 8}, {1, 6, 2}},
  {{2, 5, 1}, {9, 7, 10}},
  {{10, 6, 9}, {0, 4, 3}},
  {{5, 11, 4}, {0, 9, 1}},
  {{1, 8, 0}, {7, 10, 6}},
  {{6, 9, 7}, {3, 11, 2}},
  {{2, 10, 3}, {4, 8, 5}}
};

int facets (coord n, double alpha, coord v[12], double h)
{
  coord a[12];
  int orient[12];

  for (int i = 0; i < 12; i++) {
    coord e, d;
    double den = 0., t = alpha;
     {
      d.x = h*(cube_edge[i][0].x - 0.5);
      e.x = h*(cube_edge[i][1].x - 0.5);
      den += n.x*(e.x - d.x);
      t -= n.x*d.x;
    } 
#line 342
{
      d.y = h*(cube_edge[i][0].y - 0.5);
      e.y = h*(cube_edge[i][1].y - 0.5);
      den += n.y*(e.y - d.y);
      t -= n.y*d.y;
    } 
#line 342
{
      d.z = h*(cube_edge[i][0].z - 0.5);
      e.z = h*(cube_edge[i][1].z - 0.5);
      den += n.z*(e.z - d.z);
      t -= n.z*d.z;
    }
    orient[i] = -1;
    if (fabs (den) > 1e-10) {
      t /= den;
      if (t >= 0. && t < 1.) {
 double s = - alpha;
  {
   a[i].x = d.x + t*(e.x - d.x);
   s += n.x*e.x;
 } 
#line 353
{
   a[i].y = d.y + t*(e.y - d.y);
   s += n.y*e.y;
 } 
#line 353
{
   a[i].z = d.z + t*(e.z - d.z);
   s += n.z*e.z;
 }
 orient[i] = (s > 0.);
      }
    }
  }

  for (int i = 0; i < 12; i++) {
    int nv = 0, e = i;
    while (orient[e] >= 0) {
      int m = 0, * ne = cube_connect[e][orient[e]];
      v[nv++] = a[e];
      orient[e] = -1;
      while (m < 3 && orient[e] < 0)
 e = ne[m++];
    }
    if (nv > 2)
      return nv;
  }
  return 0;
}






double line_length_center (coord m, double alpha, coord * p)
{
  alpha += (m.x + m.y)/2.;

  coord n = m;
  
    if (n.x < 0.) {
      alpha -= n.x;
      n.x = - n.x;
    }
    
#line 388
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
    
#line 399
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
#line 424
{
    p->y /= 2.;
    p->y = clamp (p->y, 0., 1.);
    if (m.y < 0.)
      p->y = 1. - p->y;
    p->y -= 0.5;
  }

  return sqrt (ax*ax + ay*ay);
}




double plane_area_center (coord m, double alpha, coord * p)
{
  
    if (fabs (m.x) < 1e-4) {
      coord n, q;
      ((double *)&n)[0] = m.y;
      ((double *)&n)[1] = m.z;
      double length = line_length_center (n, alpha, &q);
      p->x = 0.;
      p->y = ((double *)&q)[0];
      p->z = ((double *)&q)[1];
      return length;
    }
    
#line 441
if (fabs (m.y) < 1e-4) {
      coord n, q;
      ((double *)&n)[0] = m.z;
      ((double *)&n)[1] = m.x;
      double length = line_length_center (n, alpha, &q);
      p->y = 0.;
      p->z = ((double *)&q)[0];
      p->x = ((double *)&q)[1];
      return length;
    }
    
#line 441
if (fabs (m.z) < 1e-4) {
      coord n, q;
      ((double *)&n)[0] = m.x;
      ((double *)&n)[1] = m.y;
      double length = line_length_center (n, alpha, &q);
      p->z = 0.;
      p->x = ((double *)&q)[0];
      p->y = ((double *)&q)[1];
      return length;
    }

  alpha += (m.x + m.y + m.z)/2.;
  coord n = m;
  
    if (n.x < 0.) {
      alpha -= n.x;
      n.x = - n.x;
    }
    
#line 455
if (n.y < 0.) {
      alpha -= n.y;
      n.y = - n.y;
    }
    
#line 455
if (n.z < 0.) {
      alpha -= n.z;
      n.z = - n.z;
    }

  double amax = n.x + n.y + n.z;
  if (alpha < 0. || alpha > amax) {
    p->x = p->y = p->z = 0.;
    return 0.;
  }

  double area = sq(alpha);
  p->x = p->y = p->z = area*alpha;

   {
    double b = alpha - n.x;
    if (b > 0.) {
      area -= b*b;
      p->x -= b*b*(2.*n.x + alpha);
      p->y -= b*b*b;
      p->z -= b*b*b;
    }
  } 
#line 469
{
    double b = alpha - n.y;
    if (b > 0.) {
      area -= b*b;
      p->y -= b*b*(2.*n.y + alpha);
      p->z -= b*b*b;
      p->x -= b*b*b;
    }
  } 
#line 469
{
    double b = alpha - n.z;
    if (b > 0.) {
      area -= b*b;
      p->z -= b*b*(2.*n.z + alpha);
      p->x -= b*b*b;
      p->y -= b*b*b;
    }
  }

  amax = alpha - amax;
   {
    double b = amax + n.x;
    if (b > 0.) {
      area += b*b;
      p->y += b*b*(2.*n.y + alpha - n.z);
      p->z += b*b*(2.*n.z + alpha - n.y);
      p->x += b*b*b;
    }
  } 
#line 480
{
    double b = amax + n.y;
    if (b > 0.) {
      area += b*b;
      p->z += b*b*(2.*n.z + alpha - n.x);
      p->x += b*b*(2.*n.x + alpha - n.z);
      p->y += b*b*b;
    }
  } 
#line 480
{
    double b = amax + n.z;
    if (b > 0.) {
      area += b*b;
      p->x += b*b*(2.*n.x + alpha - n.y);
      p->y += b*b*(2.*n.y + alpha - n.x);
      p->z += b*b*b;
    }
  }

  area *= 3.;
   {
    if (area) {
      p->x /= area*n.x;
      p->x = clamp (p->x, 0., 1.);
    }
    else
      p->x = 0.;
    if (m.x < 0.) p->x = 1. - p->x;
    p->x -= 0.5;
  } 
#line 491
{
    if (area) {
      p->y /= area*n.y;
      p->y = clamp (p->y, 0., 1.);
    }
    else
      p->y = 0.;
    if (m.y < 0.) p->y = 1. - p->y;
    p->y -= 0.5;
  } 
#line 491
{
    if (area) {
      p->z /= area*n.z;
      p->z = clamp (p->z, 0., 1.);
    }
    else
      p->z = 0.;
    if (m.z < 0.) p->z = 1. - p->z;
    p->z -= 0.5;
  }

  return area*sqrt (1./(sq(n.x)*sq(n.y)) +
      1./(sq(n.x)*sq(n.z)) +
      1./(sq(n.z)*sq(n.y)))/6.;
}






void line_center (coord m, double alpha, double a, coord * p)
{
  alpha += (m.x + m.y)/2.;

  coord n = m;
  
    if (n.x < 0.) {
      alpha -= n.x;
      n.x = - n.x;
    }
    
#line 518
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
    
#line 535
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
#line 543
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
#line 551
{
    p->y /= 6.*sq(n.y)*n.x*a;
    p->y = sign(m.y)*(p->y - 0.5);
  }
}
#line 564 "/home/lisergey/basilisk/src/geometry.h"
void plane_center (coord m, double alpha, double a, coord * p)
{
  
    if (fabs (m.x) < 1e-4) {
      coord n, q;
      ((double *)&n)[0] = m.y;
      ((double *)&n)[1] = m.z;
      line_center (n, alpha, a, &q);
      p->x = 0.;
      p->y = ((double *)&q)[0];
      p->z = ((double *)&q)[1];
      return;
    }
    
#line 567
if (fabs (m.y) < 1e-4) {
      coord n, q;
      ((double *)&n)[0] = m.z;
      ((double *)&n)[1] = m.x;
      line_center (n, alpha, a, &q);
      p->y = 0.;
      p->z = ((double *)&q)[0];
      p->x = ((double *)&q)[1];
      return;
    }
    
#line 567
if (fabs (m.z) < 1e-4) {
      coord n, q;
      ((double *)&n)[0] = m.x;
      ((double *)&n)[1] = m.y;
      line_center (n, alpha, a, &q);
      p->z = 0.;
      p->x = ((double *)&q)[0];
      p->y = ((double *)&q)[1];
      return;
    }

  alpha += (m.x + m.y + m.z)/2.;
  coord n = m;
  
    if (n.x < 0.) {
      alpha -= n.x;
      n.x = - n.x;
    }
    
#line 581
if (n.y < 0.) {
      alpha -= n.y;
      n.y = - n.y;
    }
    
#line 581
if (n.z < 0.) {
      alpha -= n.z;
      n.z = - n.z;
    }

  if (alpha <= 0. || a == 0.) {
    p->x = p->y = p->z = -0.5;
    return;
  }

  if (alpha >= n.x + n.y + n.z || a == 1.) {
    p->x = p->y = p->z = 0.;
    return;
  }

  p->x = p->y = p->z = sq(sq(alpha));
   {
    double b = alpha - n.x;
    if (b > 0.) {
      p->x -= cube(b)*(3.*n.x + alpha);
      p->y -= sq(sq(b));
      p->z -= sq(sq(b));
    }
  } 
#line 597
{
    double b = alpha - n.y;
    if (b > 0.) {
      p->y -= cube(b)*(3.*n.y + alpha);
      p->z -= sq(sq(b));
      p->x -= sq(sq(b));
    }
  } 
#line 597
{
    double b = alpha - n.z;
    if (b > 0.) {
      p->z -= cube(b)*(3.*n.z + alpha);
      p->x -= sq(sq(b));
      p->y -= sq(sq(b));
    }
  }

  double amax = alpha - (n.x + n.y + n.z);
   {
    double b = amax + n.z;
    if (b > 0.) {
      p->x += cube(b)*(3.*n.x + alpha - n.y);
      p->y += cube(b)*(3.*n.y + alpha - n.x);
      p->z += sq(sq(b));
    }
  } 
#line 607
{
    double b = amax + n.x;
    if (b > 0.) {
      p->y += cube(b)*(3.*n.y + alpha - n.z);
      p->z += cube(b)*(3.*n.z + alpha - n.y);
      p->x += sq(sq(b));
    }
  } 
#line 607
{
    double b = amax + n.y;
    if (b > 0.) {
      p->z += cube(b)*(3.*n.z + alpha - n.x);
      p->x += cube(b)*(3.*n.x + alpha - n.z);
      p->y += sq(sq(b));
    }
  }

  double b = 24.*n.x*n.y*n.z*a;
   {
    p->x /= b*n.x;
    p->x = sign(m.x)*(p->x - 0.5);
  } 
#line 617
{
    p->y /= b*n.y;
    p->y = sign(m.y)*(p->y - 0.5);
  } 
#line 617
{
    p->z /= b*n.z;
    p->z = sign(m.z)*(p->z - 0.5);
  }
}
#line 13 "/home/lisergey/basilisk/src/fractions.h"







#line 1 "myc.h"
#line 1 "/home/lisergey/basilisk/src/myc.h"
#line 16 "/home/lisergey/basilisk/src/myc.h"
coord mycs (Point point, scalar c)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;
  double m1,m2,m[4][3],t0,t1,t2;
  int cn;



  m1 = val(c,-1,0,-1) + val(c,-1,0,1) + val(c,-1,-1,0) + val(c,-1,1,0) +
       val(c,-1,0,0);
  m2 = val(c,1,0,-1) + val(c,1,0,1) + val(c,1,-1,0) + val(c,1,1,0) +
       val(c,1,0,0);
  m[0][0] = m1 > m2 ? 1. : -1.;

  m1 = val(c,-1,-1,0)+ val(c,1,-1,0)+ val(c,0,-1,0);
  m2 = val(c,-1,1,0)+ val(c,1,1,0)+ val(c,0,1,0);
  m[0][1] = 0.5*(m1-m2);

  m1 = val(c,-1,0,-1)+ val(c,1,0,-1)+ val(c,0,0,-1);
  m2 = val(c,-1,0,1)+ val(c,1,0,1)+ val(c,0,0,1);
  m[0][2] = 0.5*(m1-m2);



  m1 = val(c,-1,-1,0) + val(c,-1,1,0) + val(c,-1,0,0);
  m2 = val(c,1,-1,0) + val(c,1,1,0) + val(c,1,0,0);
  m[1][0] = 0.5*(m1-m2);

  m1 = val(c,0,-1,-1) + val(c,0,-1,1) + val(c,1,-1,0) + val(c,-1,-1,0) +
       val(c,0,-1,0);
  m2 = val(c,0,1,-1) + val(c,0,1,1) + val(c,1,1,0) + val(c,-1,1,0) +
       val(c,0,1,0);
  m[1][1] = m1 > m2 ? 1. : -1.;

  m1 = val(c,0,-1,-1)+ val(c,0,0,-1)+ val(c,0,1,-1);
  m2 = val(c,0,-1,1)+ val(c,0,0,1)+ val(c,0,1,1);
  m[1][2] = 0.5*(m1-m2);




  m1 = val(c,-1,0,-1)+ val(c,-1,0,1)+ val(c,-1,0,0);
  m2 = val(c,1,0,-1)+ val(c,1,0,1)+ val(c,1,0,0);
  m[2][0] = 0.5*(m1-m2);

  m1 = val(c,0,-1,-1)+ val(c,0,-1,1)+ val(c,0,-1,0);
  m2 = val(c,0,1,-1)+ val(c,0,1,1)+ val(c,0,1,0);
  m[2][1] = 0.5*(m1-m2);

  m1 = val(c,-1,0,-1) + val(c,1,0,-1) + val(c,0,-1,-1) + val(c,0,1,-1) +
       val(c,0,0,-1);
  m2 = val(c,-1,0,1) + val(c,1,0,1) + val(c,0,-1,1) + val(c,0,1,1) +
       val(c,0,0,1);
  m[2][2] = m1 > m2 ? 1. : -1.;


  t0 = fabs(m[0][0]) + fabs(m[0][1]) + fabs(m[0][2]);
  m[0][0] /= t0;
  m[0][1] /= t0;
  m[0][2] /= t0;

  t0 = fabs(m[1][0]) + fabs(m[1][1]) + fabs(m[1][2]);
  m[1][0] /= t0;
  m[1][1] /= t0;
  m[1][2] /= t0;

  t0 = fabs(m[2][0]) + fabs(m[2][1]) + fabs(m[2][2]);
  m[2][0] /= t0;
  m[2][1] /= t0;
  m[2][2] /= t0;


  t0 = fabs(m[0][0]);
  t1 = fabs(m[1][1]);
  t2 = fabs(m[2][2]);

  cn = 0;
  if (t1 > t0) {
    t0 = t1;
    cn = 1;
  }
  if (t2 > t0)
    cn = 2;


  m1 = val(c,-1,-1,-1) + val(c,-1,1,-1) + val(c,-1,-1,1) + val(c,-1,1,1) +
       2.*(val(c,-1,-1,0) + val(c,-1,1,0) + val(c,-1,0,-1) + val(c,-1,0,1)) +
       4.*val(c,-1,0,0);
  m2 = val(c,1,-1,-1) + val(c,1,1,-1) + val(c,1,-1,1) + val(c,1,1,1) +
       2.*(val(c,1,-1,0) + val(c,1,1,0) + val(c,1,0,-1) + val(c,1,0,1)) +
       4.*val(c,1,0,0);
  m[3][0] = m1 - m2;

  m1 = val(c,-1,-1,-1) + val(c,-1,-1,1) + val(c,1,-1,-1) + val(c,1,-1,1) +
       2.*( val(c,-1,-1,0) + val(c,1,-1,0) + val(c,0,-1,-1) + val(c,0,-1,1)) +
       4.*val(c,0,-1,0);
  m2 = val(c,-1,1,-1) + val(c,-1,1,1) + val(c,1,1,-1) + val(c,1,1,1) +
       2.*(val(c,-1,1,0) + val(c,1,1,0) + val(c,0,1,-1) + val(c,0,1,1)) +
       4.*val(c,0,1,0);
  m[3][1] = m1 - m2;

  m1 = val(c,-1,-1,-1) + val(c,-1,1,-1) + val(c,1,-1,-1) + val(c,1,1,-1) +
       2.*(val(c,-1,0,-1) + val(c,1,0,-1) + val(c,0,-1,-1) + val(c,0,1,-1)) +
       4.*val(c,0,0,-1);
  m2 = val(c,-1,-1,1) + val(c,-1,1,1) + val(c,1,-1,1) + val(c,1,1,1) +
       2.*(val(c,-1,0,1) + val(c,1,0,1) + val(c,0,-1,1) + val(c,0,1,1)) +
       4.*val(c,0,0,1);
  m[3][2] = m1 - m2;


  t0 = fabs(m[3][0]) + fabs(m[3][1]) + fabs(m[3][2]);
  if (t0 < 1e-30) {
    coord mxyz = {1., 0., 0.};
    return mxyz;
  }

  m[3][0] /= t0;
  m[3][1] /= t0;
  m[3][2] /= t0;


  t0 = fabs (m[3][0]);
  t1 = fabs (m[3][1]);
  t2 = fabs (m[3][2]);
  if (t1 > t0)
    t0 = t1;
  if (t2 > t0)
    t0 = t2;

  if (fabs(m[cn][cn]) > t0)
    cn = 3;


  coord mxyz = {m[cn][0], m[cn][1], m[cn][2]};
  return mxyz;
}
#line 13 "/home/lisergey/basilisk/src/fractions.h"







#line 1 "myc.h"
#line 1 "/home/lisergey/basilisk/src/myc.h"
#line 16 "/home/lisergey/basilisk/src/myc.h"
static void _stencil_mycs (Point point, scalar c)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;  
  
   



_stencil_val(c,-1,0,-1); _stencil_val(c,-1,0,1); _stencil_val(c,-1,-1,0); _stencil_val(c,-1,1,0);
       _stencil_val(c,-1,0,0); 



       
#line 25
_stencil_val(c,1,0,-1); _stencil_val(c,1,0,1); _stencil_val(c,1,-1,0); _stencil_val(c,1,1,0);
       _stencil_val(c,1,0,0);  
      
          
_stencil_val(c,-1,-1,0); _stencil_val(c,1,-1,0); _stencil_val(c,0,-1,0);

    
#line 30
_stencil_val(c,-1,1,0); _stencil_val(c,1,1,0); _stencil_val(c,0,1,0);
   
     
_stencil_val(c,-1,0,-1); _stencil_val(c,1,0,-1); _stencil_val(c,0,0,-1);

    
#line 34
_stencil_val(c,-1,0,1); _stencil_val(c,1,0,1); _stencil_val(c,0,0,1);
   
     


_stencil_val(c,-1,-1,0); _stencil_val(c,-1,1,0); _stencil_val(c,-1,0,0);



      
#line 40
_stencil_val(c,1,-1,0); _stencil_val(c,1,1,0); _stencil_val(c,1,0,0);
     
     
_stencil_val(c,0,-1,-1); _stencil_val(c,0,-1,1); _stencil_val(c,1,-1,0); _stencil_val(c,-1,-1,0);
       _stencil_val(c,0,-1,0);

        
#line 45
_stencil_val(c,0,1,-1); _stencil_val(c,0,1,1); _stencil_val(c,1,1,0); _stencil_val(c,-1,1,0);
       _stencil_val(c,0,1,0);
       
           
_stencil_val(c,0,-1,-1); _stencil_val(c,0,0,-1); _stencil_val(c,0,1,-1);

    
#line 50
_stencil_val(c,0,-1,1); _stencil_val(c,0,0,1); _stencil_val(c,0,1,1);
   
     



_stencil_val(c,-1,0,-1); _stencil_val(c,-1,0,1); _stencil_val(c,-1,0,0);




    
#line 57
_stencil_val(c,1,0,-1); _stencil_val(c,1,0,1); _stencil_val(c,1,0,0);
   
     
_stencil_val(c,0,-1,-1); _stencil_val(c,0,-1,1); _stencil_val(c,0,-1,0);

    
#line 61
_stencil_val(c,0,1,-1); _stencil_val(c,0,1,1); _stencil_val(c,0,1,0);
   
     
_stencil_val(c,-1,0,-1); _stencil_val(c,1,0,-1); _stencil_val(c,0,-1,-1); _stencil_val(c,0,1,-1);
       _stencil_val(c,0,0,-1);

        
#line 66
_stencil_val(c,-1,0,1); _stencil_val(c,1,0,1); _stencil_val(c,0,-1,1); _stencil_val(c,0,1,1);
       _stencil_val(c,0,0,1);    
       
          


       
    
    
    

        
    
    
    

        
    
    
    


    
   
   

    
      
     
   
      
     
      

_stencil_val(c,-1,-1,-1); _stencil_val(c,-1,1,-1); _stencil_val(c,-1,-1,1); _stencil_val(c,-1,1,1);
_stencil_val(c,-1,-1,0); _stencil_val(c,-1,1,0); _stencil_val(c,-1,0,-1); _stencil_val(c,-1,0,1);
_stencil_val(c,-1,0,0);


       
           
        
#line 103
_stencil_val(c,1,-1,-1); _stencil_val(c,1,1,-1); _stencil_val(c,1,-1,1); _stencil_val(c,1,1,1);
_stencil_val(c,1,-1,0); _stencil_val(c,1,1,0); _stencil_val(c,1,0,-1); _stencil_val(c,1,0,1);
_stencil_val(c,1,0,0);
       
           
       
       
#line 108
_stencil_val(c,-1,-1,-1); _stencil_val(c,-1,-1,1); _stencil_val(c,1,-1,-1); _stencil_val(c,1,-1,1); 
_stencil_val(c,-1,-1,0); _stencil_val(c,1,-1,0); _stencil_val(c,0,-1,-1); _stencil_val(c,0,-1,1);
_stencil_val(c,0,-1,0);

       
           
        
#line 111
_stencil_val(c,-1,1,-1); _stencil_val(c,-1,1,1); _stencil_val(c,1,1,-1); _stencil_val(c,1,1,1);
_stencil_val(c,-1,1,0); _stencil_val(c,1,1,0); _stencil_val(c,0,1,-1); _stencil_val(c,0,1,1);
_stencil_val(c,0,1,0);
       
           
       
       
#line 116
_stencil_val(c,-1,-1,-1); _stencil_val(c,-1,1,-1); _stencil_val(c,1,-1,-1); _stencil_val(c,1,1,-1);
_stencil_val(c,-1,0,-1); _stencil_val(c,1,0,-1); _stencil_val(c,0,-1,-1); _stencil_val(c,0,1,-1);
_stencil_val(c,0,0,-1);

       
           
        
#line 119
_stencil_val(c,-1,-1,1); _stencil_val(c,-1,1,1); _stencil_val(c,1,-1,1); _stencil_val(c,1,1,1);
_stencil_val(c,-1,0,1); _stencil_val(c,1,0,1); _stencil_val(c,0,-1,1); _stencil_val(c,0,1,1);
_stencil_val(c,0,0,1);          
     
    
          
       
           
       
      


        
     

    
    
    


     
     
     
  
      
  
      
     

     


  
  
#line 149
return ;
}
#line 21 "/home/lisergey/basilisk/src/fractions.h"
#line 40 "/home/lisergey/basilisk/src/fractions.h"
void fraction_refine (Point point, scalar c)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;





  double cc = val(c,0,0,0);
  if (cc <= 0. || cc >= 1.)
    {foreach_child()
      val(c,0,0,0) = cc;end_foreach_child()}
  else {




    coord n = mycs (point, c);
    double alpha = plane_alpha (cc, n);






    {foreach_child() {
      static const coord a = {0.,0.,0.}, b = {.5,.5,.5};
      coord nc;
      
 nc.x = child.x*n.x;
 
#line 68
nc.y = child.y*n.y;
 
#line 68
nc.z = child.z*n.z;
      val(c,0,0,0) = rectangle_fraction (nc, alpha, a, b);
    }end_foreach_child()}
  }
}











static void alpha_refine (Point point, scalar alpha)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;
  vector n = _attribute[alpha.i].n;
  double alphac = 2.*val(alpha,0,0,0);
  coord m;
  
    m.x = val(n.x,0,0,0);
    
#line 90
m.y = val(n.y,0,0,0);
    
#line 90
m.z = val(n.z,0,0,0);
  {foreach_child() {
    val(alpha,0,0,0) = alphac;
    
      val(alpha,0,0,0) -= child.x*m.x/2.;
      
#line 94
val(alpha,0,0,0) -= child.y*m.y/2.;
      
#line 94
val(alpha,0,0,0) -= child.z*m.z/2.;
  }end_foreach_child()}
}
#line 120 "/home/lisergey/basilisk/src/fractions.h"
     
void fractions (scalar Phi, scalar c,
  vector s, double val)
{tracing("fractions","/home/lisergey/basilisk/src/fractions.h",121);

  vector   as=(s).x.i>0?(s):new_face_vector("as");
#line 134 "/home/lisergey/basilisk/src/fractions.h"
  vector  p=new_vector("p");
#line 146 "/home/lisergey/basilisk/src/fractions.h"
  foreach_vertex_stencil(1,{(NonLocal[]){{"p","vector",(void *)&p,NULL,0},{"Phi","scalar",(void *)&Phi,NULL,0},{"val","double",(void *)&val,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},""," \n// #line 146 \"/home/lisergey/basilisk/src/fractions.h\"\n{ if (is_vertex(neighbor(1,0,0))) {\n\n\n\n\n\n    if ((val(Phi,0,0,0) - val)*(val(Phi,1,0,0) - val) < 0.) {\n\n\n\n\n\n\n      val_out_(p.x,0,0,0) = (val(Phi,0,0,0) - val)/(val(Phi,0,0,0) - val(Phi,1,0,0));\n      if (val(Phi,0,0,0) < val)\n val_out_(p.x,0,0,0) = 1. - val(p.x,0,0,0);\n    }\n// # 171 \"/home/lisergey/basilisk/src/fractions.h\"\n    else\n      val_out_(p.x,0,0,0) = (val(Phi,0,0,0) > val || val(Phi,1,0,0) > val);\n  } \n// #line 146\nif (is_vertex(neighbor(0,1,0))) {\n\n\n\n\n\n    if ((val(Phi,0,0,0) - val)*(val(Phi,0,1,0) - val) < 0.) {\n\n\n\n\n\n\n      val_out_(p.y,0,0,0) = (val(Phi,0,0,0) - val)/(val(Phi,0,0,0) - val(Phi,0,1,0));\n      if (val(Phi,0,0,0) < val)\n val_out_(p.y,0,0,0) = 1. - val(p.y,0,0,0);\n    }\n// # 171 \"/home/lisergey/basilisk/src/fractions.h\"\n    else\n      val_out_(p.y,0,0,0) = (val(Phi,0,0,0) > val || val(Phi,0,1,0) > val);\n  } \n// #line 146\nif (is_vertex(neighbor(0,0,1))) {\n\n\n\n\n\n    if ((val(Phi,0,0,0) - val)*(val(Phi,0,0,1) - val) < 0.) {\n\n\n\n\n\n\n      val_out_(p.z,0,0,0) = (val(Phi,0,0,0) - val)/(val(Phi,0,0,0) - val(Phi,0,0,1));\n      if (val(Phi,0,0,0) < val)\n val_out_(p.z,0,0,0) = 1. - val(p.z,0,0,0);\n    }\n// # 171 \"/home/lisergey/basilisk/src/fractions.h\"\n    else\n      val_out_(p.z,0,0,0) = (val(Phi,0,0,0) > val || val(Phi,0,0,1) > val);\n  }}"}) { {_stencil_neighbor(1,0,0); {





_stencil_val(Phi,0,0,0);_stencil_val(Phi,1,0,0);{ {






_stencil_val(Phi,0,0,0);_stencil_val(Phi,0,0,0); _stencil_val(Phi,1,0,0);






      
#line 159
_stencil_val_a(p.x,0,0,0);
_stencil_val(Phi,0,0,0);
 { _stencil_val(p.x,0,0,0);_stencil_val_a(p.x,0,0,0);   }     
         
    
#line 162
}
      








{_stencil_val(Phi,0,0,0); _stencil_val(Phi,1,0,0);_stencil_val_a(p.x,0,0,0);       }}





           
#line 171 "/home/lisergey/basilisk/src/fractions.h"
    
  
} } 
#line 146
{_stencil_neighbor(0,1,0); {





_stencil_val(Phi,0,0,0);_stencil_val(Phi,0,1,0);{ {






_stencil_val(Phi,0,0,0);_stencil_val(Phi,0,0,0); _stencil_val(Phi,0,1,0);






      
#line 159
_stencil_val_a(p.y,0,0,0);
_stencil_val(Phi,0,0,0);
 { _stencil_val(p.y,0,0,0);_stencil_val_a(p.y,0,0,0);   }     
         
    
#line 162
}
      








{_stencil_val(Phi,0,0,0); _stencil_val(Phi,0,1,0);_stencil_val_a(p.y,0,0,0);       }}





           
#line 171 "/home/lisergey/basilisk/src/fractions.h"
    
  
} } 
#line 146
{_stencil_neighbor(0,0,1); {





_stencil_val(Phi,0,0,0);_stencil_val(Phi,0,0,1);{ {






_stencil_val(Phi,0,0,0);_stencil_val(Phi,0,0,0); _stencil_val(Phi,0,0,1);






      
#line 159
_stencil_val_a(p.z,0,0,0);
_stencil_val(Phi,0,0,0);
 { _stencil_val(p.z,0,0,0);_stencil_val_a(p.z,0,0,0);   }     
         
    
#line 162
}
      








{_stencil_val(Phi,0,0,0); _stencil_val(Phi,0,0,1);_stencil_val_a(p.z,0,0,0);       }}





           
#line 171 "/home/lisergey/basilisk/src/fractions.h"
    
  
} }}end_foreach_vertex_stencil();
#line 146 "/home/lisergey/basilisk/src/fractions.h"
  {foreach_vertex() { if (is_vertex(neighbor(1,0,0))) {





    if ((val(Phi,0,0,0) - val)*(val(Phi,1,0,0) - val) < 0.) {






      val(p.x,0,0,0) = (val(Phi,0,0,0) - val)/(val(Phi,0,0,0) - val(Phi,1,0,0));
      if (val(Phi,0,0,0) < val)
 val(p.x,0,0,0) = 1. - val(p.x,0,0,0);
    }
#line 171 "/home/lisergey/basilisk/src/fractions.h"
    else
      val(p.x,0,0,0) = (val(Phi,0,0,0) > val || val(Phi,1,0,0) > val);
  } 
#line 146
if (is_vertex(neighbor(0,1,0))) {





    if ((val(Phi,0,0,0) - val)*(val(Phi,0,1,0) - val) < 0.) {






      val(p.y,0,0,0) = (val(Phi,0,0,0) - val)/(val(Phi,0,0,0) - val(Phi,0,1,0));
      if (val(Phi,0,0,0) < val)
 val(p.y,0,0,0) = 1. - val(p.y,0,0,0);
    }
#line 171 "/home/lisergey/basilisk/src/fractions.h"
    else
      val(p.y,0,0,0) = (val(Phi,0,0,0) > val || val(Phi,0,1,0) > val);
  } 
#line 146
if (is_vertex(neighbor(0,0,1))) {





    if ((val(Phi,0,0,0) - val)*(val(Phi,0,0,1) - val) < 0.) {






      val(p.z,0,0,0) = (val(Phi,0,0,0) - val)/(val(Phi,0,0,0) - val(Phi,0,0,1));
      if (val(Phi,0,0,0) < val)
 val(p.z,0,0,0) = 1. - val(p.z,0,0,0);
    }
#line 171 "/home/lisergey/basilisk/src/fractions.h"
    else
      val(p.z,0,0,0) = (val(Phi,0,0,0) > val || val(Phi,0,0,1) > val);
  }}end_foreach_vertex();}
#line 190 "/home/lisergey/basilisk/src/fractions.h"
  
    _attribute[p.x.i].dirty = false;
    
#line 191
_attribute[p.y.i].dirty = false;
    
#line 191
_attribute[p.z.i].dirty = false;

  scalar s_x = as.x, s_y = as.y, s_z = as.z;
  foreach_face_stencil(1,{(NonLocal[]){{"s_y","scalar",(void *)&s_y,NULL,0},{"s_x","scalar",(void *)&s_x,NULL,0},{"Phi","scalar",(void *)&Phi,NULL,0},{"s_z","scalar",(void *)&s_z,NULL,0},{"p","vector",(void *)&p,NULL,0},{"val","double",(void *)&val,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"\n// # 163 \"/home/lisergey/basilisk/src/geometry.h\"\nreal line_area (real nx, real ny, real alpha)\n{\n  real a, v, area;\n\n  alpha += (nx + ny)/2.;\n  if (nx < 0.) {\n    alpha -= nx;\n    nx = - nx;\n  }\n  if (ny < 0.) {\n    alpha -= ny;\n    ny = - ny;\n  }\n\n  if (alpha <= 0.)\n    return 0.;\n\n  if (alpha >= nx + ny)\n    return 1.;\n\n  if (nx < 1e-10)\n    area = alpha/ny;\n  else if (ny < 1e-10)\n    area = alpha/nx;\n  else {\n    v = sq(alpha);\n\n    a = alpha - nx;\n    if (a > 0.)\n      v -= a*a;\n\n    a = alpha - ny;\n    if (a > 0.)\n      v -= a*a;\n\n    area = v/(2.*nx*ny);\n  }\n\n  return clamp (area, 0., 1.);\n}","\n// #line 194 \"/home/lisergey/basilisk/src/fractions.h\"\n{is_face_z(){\n\n\n\n\n  {\n// # 231 \"/home/lisergey/basilisk/src/fractions.h\"\n    coord n;\n    real nn = 0.;\n     {\n      n.x = val(p.y,0,0,0) - val(p.y,1,0,0);\n      nn += fabs(n.x);\n    } \n// #line 233\n{\n      n.y = val(p.x,0,0,0) - val(p.x,0,1,0);\n      nn += fabs(n.y);\n    }\n\n\n\n\n\n    if (nn == 0.)\n      val_out_(s_z,0,0,0) = val(p.x,0,0,0);\n    else {\n\n\n\n\n\n      \n n.x /= nn;\n \n// #line 251\nn.y /= nn;\n\n\n\n\n\n\n      real alpha = 0., ni = 0.;\n      for (int i = 0; i <= 1; i++)\n {\n   if (val(p.x,0,i,0) > 0. && val(p.x,0,i,0) < 1.) {\n     real a = sign(val(Phi,0,i,0) - val)*(val(p.x,0,i,0) - 0.5);\n     alpha += n.x*a + n.y*(i - 0.5);\n     ni++;\n   }\n   \n// #line 261\nif (val(p.y,i,0,0) > 0. && val(p.y,i,0,0) < 1.) {\n     real a = sign(val(Phi,i,0,0) - val)*(val(p.y,i,0,0) - 0.5);\n     alpha += n.y*a + n.x*(i - 0.5);\n     ni++;\n   }}\n// # 274 \"/home/lisergey/basilisk/src/fractions.h\"\n      if (ni == 0)\n val_out_(s_z,0,0,0) = max (val(p.x,0,0,0), val(p.y,0,0,0));\n      else if (ni != 4)\n val_out_(s_z,0,0,0) = line_area (n.x, n.y, alpha/ni);\n      else {\n\n val_out_(s_z,0,0,0) = (val(p.x,0,0,0) + val(p.x,0,1,0) + val(p.y,0,0,0) + val(p.y,1,0,0) > 2.);\n\n\n\n      }\n    }\n  }}end_is_face_z()\n// #line 194\nis_face_x(){\n\n\n\n\n  {\n// # 231 \"/home/lisergey/basilisk/src/fractions.h\"\n    coord n;\n    real nn = 0.;\n     {\n      n.y = val(p.z,0,0,0) - val(p.z,0,1,0);\n      nn += fabs(n.y);\n    } \n// #line 233\n{\n      n.z = val(p.y,0,0,0) - val(p.y,0,0,1);\n      nn += fabs(n.z);\n    }\n\n\n\n\n\n    if (nn == 0.)\n      val_out_(s_x,0,0,0) = val(p.y,0,0,0);\n    else {\n\n\n\n\n\n      \n n.y /= nn;\n \n// #line 251\nn.z /= nn;\n\n\n\n\n\n\n      real alpha = 0., ni = 0.;\n      for (int i = 0; i <= 1; i++)\n {\n   if (val(p.y,0,0,i) > 0. && val(p.y,0,0,i) < 1.) {\n     real a = sign(val(Phi,0,0,i) - val)*(val(p.y,0,0,i) - 0.5);\n     alpha += n.y*a + n.z*(i - 0.5);\n     ni++;\n   }\n   \n// #line 261\nif (val(p.z,0,i,0) > 0. && val(p.z,0,i,0) < 1.) {\n     real a = sign(val(Phi,0,i,0) - val)*(val(p.z,0,i,0) - 0.5);\n     alpha += n.z*a + n.y*(i - 0.5);\n     ni++;\n   }}\n// # 274 \"/home/lisergey/basilisk/src/fractions.h\"\n      if (ni == 0)\n val_out_(s_x,0,0,0) = max (val(p.y,0,0,0), val(p.z,0,0,0));\n      else if (ni != 4)\n val_out_(s_x,0,0,0) = line_area (n.y, n.z, alpha/ni);\n      else {\n\n val_out_(s_x,0,0,0) = (val(p.y,0,0,0) + val(p.y,0,0,1) + val(p.z,0,0,0) + val(p.z,0,1,0) > 2.);\n\n\n\n      }\n    }\n  }}end_is_face_x()\n// #line 194\nis_face_y(){\n\n\n\n\n  {\n// # 231 \"/home/lisergey/basilisk/src/fractions.h\"\n    coord n;\n    real nn = 0.;\n     {\n      n.z = val(p.x,0,0,0) - val(p.x,0,0,1);\n      nn += fabs(n.z);\n    } \n// #line 233\n{\n      n.x = val(p.z,0,0,0) - val(p.z,1,0,0);\n      nn += fabs(n.x);\n    }\n\n\n\n\n\n    if (nn == 0.)\n      val_out_(s_y,0,0,0) = val(p.z,0,0,0);\n    else {\n\n\n\n\n\n      \n n.z /= nn;\n \n// #line 251\nn.x /= nn;\n\n\n\n\n\n\n      real alpha = 0., ni = 0.;\n      for (int i = 0; i <= 1; i++)\n {\n   if (val(p.z,i,0,0) > 0. && val(p.z,i,0,0) < 1.) {\n     real a = sign(val(Phi,i,0,0) - val)*(val(p.z,i,0,0) - 0.5);\n     alpha += n.z*a + n.x*(i - 0.5);\n     ni++;\n   }\n   \n// #line 261\nif (val(p.x,0,0,i) > 0. && val(p.x,0,0,i) < 1.) {\n     real a = sign(val(Phi,0,0,i) - val)*(val(p.x,0,0,i) - 0.5);\n     alpha += n.x*a + n.z*(i - 0.5);\n     ni++;\n   }}\n// # 274 \"/home/lisergey/basilisk/src/fractions.h\"\n      if (ni == 0)\n val_out_(s_y,0,0,0) = max (val(p.z,0,0,0), val(p.x,0,0,0));\n      else if (ni != 4)\n val_out_(s_y,0,0,0) = line_area (n.z, n.x, alpha/ni);\n      else {\n\n val_out_(s_y,0,0,0) = (val(p.z,0,0,0) + val(p.z,1,0,0) + val(p.x,0,0,0) + val(p.x,0,0,1) > 2.);\n\n\n\n      }\n    }\n  }}end_is_face_y()}"}){_stencil_is_face_z(){




  {    
#line 231 "/home/lisergey/basilisk/src/fractions.h"
    
    
     { 
_stencil_val(p.y,0,0,0); _stencil_val(p.y,1,0,0);  
       
       
    
#line 236
} 
#line 233
{ 
_stencil_val(p.x,0,0,0); _stencil_val(p.x,0,1,0);  
       
       
    
#line 236
}





{
      { _stencil_val(p.x,0,0,0);_stencil_val_a(s_z,0,0,0); } 
{      





      
   






      
      for (int i = 0; i <= 1; i++)
 {
   {_stencil_val(p.x,0,i,0); _stencil_val(p.x,0,i,0); {       
     _stencil_val(p.x,0,i,0);_stencil_val(Phi,0,i,0); 
          
     
   }      }
   
#line 261
{_stencil_val(p.y,i,0,0); _stencil_val(p.y,i,0,0); {       
     _stencil_val(p.y,i,0,0);_stencil_val(Phi,i,0,0); 
          
     
   }      }}








{
 {_stencil_val(p.x,0,0,0); _stencil_val(p.y,0,0,0);_stencil_val_a(s_z,0,0,0);   }
{
 {_stencil_val_a(s_z,0,0,0);     } 
{

_stencil_val(p.x,0,0,0); _stencil_val(p.x,0,1,0); _stencil_val(p.y,0,0,0); _stencil_val(p.y,1,0,0);

 
#line 280
_stencil_val_a(s_z,0,0,0);       



      }}}
#line 274 "/home/lisergey/basilisk/src/fractions.h"
         
          
      
    







}}





       
    
  
#line 286
}}end__stencil_is_face_z()
#line 194
_stencil_is_face_x(){




  {    
#line 231 "/home/lisergey/basilisk/src/fractions.h"
    
    
     { 
_stencil_val(p.z,0,0,0); _stencil_val(p.z,0,1,0);  
       
       
    
#line 236
} 
#line 233
{ 
_stencil_val(p.y,0,0,0); _stencil_val(p.y,0,0,1);  
       
       
    
#line 236
}





{
      { _stencil_val(p.y,0,0,0);_stencil_val_a(s_x,0,0,0); } 
{      





      
   






      
      for (int i = 0; i <= 1; i++)
 {
   {_stencil_val(p.y,0,0,i); _stencil_val(p.y,0,0,i); {       
     _stencil_val(p.y,0,0,i);_stencil_val(Phi,0,0,i); 
          
     
   }      }
   
#line 261
{_stencil_val(p.z,0,i,0); _stencil_val(p.z,0,i,0); {       
     _stencil_val(p.z,0,i,0);_stencil_val(Phi,0,i,0); 
          
     
   }      }}








{
 {_stencil_val(p.y,0,0,0); _stencil_val(p.z,0,0,0);_stencil_val_a(s_x,0,0,0);   }
{
 {_stencil_val_a(s_x,0,0,0);     } 
{

_stencil_val(p.y,0,0,0); _stencil_val(p.y,0,0,1); _stencil_val(p.z,0,0,0); _stencil_val(p.z,0,1,0);

 
#line 280
_stencil_val_a(s_x,0,0,0);       



      }}}
#line 274 "/home/lisergey/basilisk/src/fractions.h"
         
          
      
    







}}





       
    
  
#line 286
}}end__stencil_is_face_x()
#line 194
_stencil_is_face_y(){




  {    
#line 231 "/home/lisergey/basilisk/src/fractions.h"
    
    
     { 
_stencil_val(p.x,0,0,0); _stencil_val(p.x,0,0,1);  
       
       
    
#line 236
} 
#line 233
{ 
_stencil_val(p.z,0,0,0); _stencil_val(p.z,1,0,0);  
       
       
    
#line 236
}





{
      { _stencil_val(p.z,0,0,0);_stencil_val_a(s_y,0,0,0); } 
{      





      
   






      
      for (int i = 0; i <= 1; i++)
 {
   {_stencil_val(p.z,i,0,0); _stencil_val(p.z,i,0,0); {       
     _stencil_val(p.z,i,0,0);_stencil_val(Phi,i,0,0); 
          
     
   }      }
   
#line 261
{_stencil_val(p.x,0,0,i); _stencil_val(p.x,0,0,i); {       
     _stencil_val(p.x,0,0,i);_stencil_val(Phi,0,0,i); 
          
     
   }      }}








{
 {_stencil_val(p.z,0,0,0); _stencil_val(p.x,0,0,0);_stencil_val_a(s_y,0,0,0);   }
{
 {_stencil_val_a(s_y,0,0,0);     } 
{

_stencil_val(p.z,0,0,0); _stencil_val(p.z,1,0,0); _stencil_val(p.x,0,0,0); _stencil_val(p.x,0,0,1);

 
#line 280
_stencil_val_a(s_y,0,0,0);       



      }}}
#line 274 "/home/lisergey/basilisk/src/fractions.h"
         
          
      
    







}}





       
    
  
#line 286
}}end__stencil_is_face_y()}end_foreach_face_stencil();
  {
#line 194
foreach_face_generic(){is_face_z(){




  {
#line 231 "/home/lisergey/basilisk/src/fractions.h"
    coord n;
    double nn = 0.;
     {
      n.x = val(p.y,0,0,0) - val(p.y,1,0,0);
      nn += fabs(n.x);
    } 
#line 233
{
      n.y = val(p.x,0,0,0) - val(p.x,0,1,0);
      nn += fabs(n.y);
    }





    if (nn == 0.)
      val(s_z,0,0,0) = val(p.x,0,0,0);
    else {





      
 n.x /= nn;
 
#line 251
n.y /= nn;






      double alpha = 0., ni = 0.;
      for (int i = 0; i <= 1; i++)
 {
   if (val(p.x,0,i,0) > 0. && val(p.x,0,i,0) < 1.) {
     double a = sign(val(Phi,0,i,0) - val)*(val(p.x,0,i,0) - 0.5);
     alpha += n.x*a + n.y*(i - 0.5);
     ni++;
   }
   
#line 261
if (val(p.y,i,0,0) > 0. && val(p.y,i,0,0) < 1.) {
     double a = sign(val(Phi,i,0,0) - val)*(val(p.y,i,0,0) - 0.5);
     alpha += n.y*a + n.x*(i - 0.5);
     ni++;
   }}
#line 274 "/home/lisergey/basilisk/src/fractions.h"
      if (ni == 0)
 val(s_z,0,0,0) = max (val(p.x,0,0,0), val(p.y,0,0,0));
      else if (ni != 4)
 val(s_z,0,0,0) = line_area (n.x, n.y, alpha/ni);
      else {

 val(s_z,0,0,0) = (val(p.x,0,0,0) + val(p.x,0,1,0) + val(p.y,0,0,0) + val(p.y,1,0,0) > 2.);



      }
    }
  }}end_is_face_z()
#line 194
is_face_x(){




  {
#line 231 "/home/lisergey/basilisk/src/fractions.h"
    coord n;
    double nn = 0.;
     {
      n.y = val(p.z,0,0,0) - val(p.z,0,1,0);
      nn += fabs(n.y);
    } 
#line 233
{
      n.z = val(p.y,0,0,0) - val(p.y,0,0,1);
      nn += fabs(n.z);
    }





    if (nn == 0.)
      val(s_x,0,0,0) = val(p.y,0,0,0);
    else {





      
 n.y /= nn;
 
#line 251
n.z /= nn;






      double alpha = 0., ni = 0.;
      for (int i = 0; i <= 1; i++)
 {
   if (val(p.y,0,0,i) > 0. && val(p.y,0,0,i) < 1.) {
     double a = sign(val(Phi,0,0,i) - val)*(val(p.y,0,0,i) - 0.5);
     alpha += n.y*a + n.z*(i - 0.5);
     ni++;
   }
   
#line 261
if (val(p.z,0,i,0) > 0. && val(p.z,0,i,0) < 1.) {
     double a = sign(val(Phi,0,i,0) - val)*(val(p.z,0,i,0) - 0.5);
     alpha += n.z*a + n.y*(i - 0.5);
     ni++;
   }}
#line 274 "/home/lisergey/basilisk/src/fractions.h"
      if (ni == 0)
 val(s_x,0,0,0) = max (val(p.y,0,0,0), val(p.z,0,0,0));
      else if (ni != 4)
 val(s_x,0,0,0) = line_area (n.y, n.z, alpha/ni);
      else {

 val(s_x,0,0,0) = (val(p.y,0,0,0) + val(p.y,0,0,1) + val(p.z,0,0,0) + val(p.z,0,1,0) > 2.);



      }
    }
  }}end_is_face_x()
#line 194
is_face_y(){




  {
#line 231 "/home/lisergey/basilisk/src/fractions.h"
    coord n;
    double nn = 0.;
     {
      n.z = val(p.x,0,0,0) - val(p.x,0,0,1);
      nn += fabs(n.z);
    } 
#line 233
{
      n.x = val(p.z,0,0,0) - val(p.z,1,0,0);
      nn += fabs(n.x);
    }





    if (nn == 0.)
      val(s_y,0,0,0) = val(p.z,0,0,0);
    else {





      
 n.z /= nn;
 
#line 251
n.x /= nn;






      double alpha = 0., ni = 0.;
      for (int i = 0; i <= 1; i++)
 {
   if (val(p.z,i,0,0) > 0. && val(p.z,i,0,0) < 1.) {
     double a = sign(val(Phi,i,0,0) - val)*(val(p.z,i,0,0) - 0.5);
     alpha += n.z*a + n.x*(i - 0.5);
     ni++;
   }
   
#line 261
if (val(p.x,0,0,i) > 0. && val(p.x,0,0,i) < 1.) {
     double a = sign(val(Phi,0,0,i) - val)*(val(p.x,0,0,i) - 0.5);
     alpha += n.x*a + n.z*(i - 0.5);
     ni++;
   }}
#line 274 "/home/lisergey/basilisk/src/fractions.h"
      if (ni == 0)
 val(s_y,0,0,0) = max (val(p.z,0,0,0), val(p.x,0,0,0));
      else if (ni != 4)
 val(s_y,0,0,0) = line_area (n.z, n.x, alpha/ni);
      else {

 val(s_y,0,0,0) = (val(p.z,0,0,0) + val(p.z,1,0,0) + val(p.x,0,0,0) + val(p.x,0,0,1) > 2.);



      }
    }
  }}end_is_face_y()}end_foreach_face_generic();}







  foreach_stencil(1,{(NonLocal[]){{"Phi","scalar",(void *)&Phi,NULL,0},{"p","vector",(void *)&p,NULL,0},{"c","scalar",(void *)&c,NULL,0},{"as","vector",(void *)&as,NULL,0},{"val","double",(void *)&val,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"\n\n\n\n\n// #line 206 \"/home/lisergey/basilisk/src/geometry.h\"\nreal plane_volume (coord n, real alpha)\n{\n  real al = alpha + (n.x + n.y + n.z)/2. +\n    max(0., -n.x) + max(0., -n.y) + max(0., -n.z);\n  if (al <= 0.)\n    return 0.;\n  real tmp = fabs(n.x) + fabs(n.y) + fabs(n.z);\n  if (al >= tmp)\n    return 1.;\n  if (tmp < 1e-10)\n    return 0.;\n  real n1 = fabs(n.x)/tmp;\n  real n2 = fabs(n.y)/tmp;\n  real n3 = fabs(n.z)/tmp;\n  al = max(0., min(1., al/tmp));\n  real al0 = min(al, 1. - al);\n  real b1 = min(n1, n2);\n  real b3 = max(n1, n2);\n  real b2 = n3;\n  if (b2 < b1) {\n    tmp = b1;\n    b1 = b2;\n    b2 = tmp;\n  }\n  else if (b2 > b3) {\n    tmp = b3;\n    b3 = b2;\n    b2 = tmp;\n  }\n  real b12 = b1 + b2;\n  real bm = min(b12, b3);\n  real pr = max(6.*b1*b2*b3, 1e-50);\n  if (al0 < b1)\n    tmp = al0*al0*al0/pr;\n  else if (al0 < b2)\n    tmp = 0.5*al0*(al0 - b1)/(b2*b3) + b1*b1*b1/pr;\n  else if (al0 < bm)\n    tmp = (al0*al0*(3.*b12 - al0) + b1*b1*(b1 - 3.*al0) +\n    b2*b2*(b2 - 3.*al0))/pr;\n  else if (b12 < b3)\n    tmp = (al0 - 0.5*bm)/b3;\n  else\n    tmp = (al0*al0*(3. - 2.*al0) + b1*b1*(b1 - 3.*al0) +\n    b2*b2*(b2 - 3.*al0) + b3*b3*(b3 - 3.*al0))/pr;\n\n  real volume = al <= 0.5 ? tmp : 1. - tmp;\n  return clamp (volume, 0., 1.);\n}"," \n// #line 294 \"/home/lisergey/basilisk/src/fractions.h\"\n{\n\n\n\n\n    coord n;\n    real nn = 0.;\n     {\n      n.x = val(as.x,0,0,0) - val(as.x,1,0,0);\n      nn += fabs(n.x);\n    } \n// #line 301\n{\n      n.y = val(as.y,0,0,0) - val(as.y,0,1,0);\n      nn += fabs(n.y);\n    } \n// #line 301\n{\n      n.z = val(as.z,0,0,0) - val(as.z,0,0,1);\n      nn += fabs(n.z);\n    }\n    if (nn == 0.)\n      val_out_(c,0,0,0) = val(as.x,0,0,0);\n    else {\n      \n n.x /= nn;\n \n// #line 309\nn.y /= nn;\n \n// #line 309\nn.z /= nn;\n\n\n\n\n\n\n      real alpha = 0., ni = 0.;\n      for (int i = 0; i <= 1; i++)\n for (int j = 0; j <= 1; j++)\n   {\n     if (val(p.x,0,i,j) > 0. && val(p.x,0,i,j) < 1.) {\n       real a = sign(val(Phi,0,i,j) - val)*(val(p.x,0,i,j) - 0.5);\n       alpha += n.x*a + n.y*(i - 0.5) + n.z*(j - 0.5);\n       ni++;\n     }\n     \n// #line 320\nif (val(p.y,j,0,i) > 0. && val(p.y,j,0,i) < 1.) {\n       real a = sign(val(Phi,j,0,i) - val)*(val(p.y,j,0,i) - 0.5);\n       alpha += n.y*a + n.z*(i - 0.5) + n.x*(j - 0.5);\n       ni++;\n     }\n     \n// #line 320\nif (val(p.z,i,j,0) > 0. && val(p.z,i,j,0) < 1.) {\n       real a = sign(val(Phi,i,j,0) - val)*(val(p.z,i,j,0) - 0.5);\n       alpha += n.z*a + n.x*(i - 0.5) + n.y*(j - 0.5);\n       ni++;\n     }}\n\n\n\n\n      if (ni == 0)\n val_out_(c,0,0,0) = val(as.x,0,0,0);\n      else if (ni < 3 || ni > 6)\n val_out_(c,0,0,0) = 0.;\n      else\n val_out_(c,0,0,0) = plane_volume (n, alpha/ni);\n    }\n  }"}) {    




    
    
     { 
_stencil_val(as.x,0,0,0); _stencil_val(as.x,1,0,0);  
       
       
    
#line 304
} 
#line 301
{ 
_stencil_val(as.y,0,0,0); _stencil_val(as.y,0,1,0);  
       
       
    
#line 304
} 
#line 301
{ 
_stencil_val(as.z,0,0,0); _stencil_val(as.z,0,0,1);  
       
       
    
#line 304
}
{
      { _stencil_val(as.x,0,0,0);_stencil_val_a(c,0,0,0); } 
{      
      
   






      
      for (int i = 0; i <= 1; i++)
 for (int j = 0; j <= 1; j++)
   {
     {_stencil_val(p.x,0,i,j); _stencil_val(p.x,0,i,j); {       
       _stencil_val(p.x,0,i,j);_stencil_val(Phi,0,i,j); 
                
       
     }      }
     
#line 320
{_stencil_val(p.y,j,0,i); _stencil_val(p.y,j,0,i); {       
       _stencil_val(p.y,j,0,i);_stencil_val(Phi,j,0,i); 
                
       
     }      }
     
#line 320
{_stencil_val(p.z,i,j,0); _stencil_val(p.z,i,j,0); {       
       _stencil_val(p.z,i,j,0);_stencil_val(Phi,i,j,0); 
                
       
     }      }}




{
 { _stencil_val(as.x,0,0,0);_stencil_val_a(c,0,0,0); }
{
 {_stencil_val_a(c,0,0,0);  }
 
{_stencil_val_a(c,0,0,0);    }}}




         
              
      
    
#line 335
}}
       
    
  
#line 336
}end_foreach_stencil();







  {
#line 294
foreach() {




    coord n;
    double nn = 0.;
     {
      n.x = val(as.x,0,0,0) - val(as.x,1,0,0);
      nn += fabs(n.x);
    } 
#line 301
{
      n.y = val(as.y,0,0,0) - val(as.y,0,1,0);
      nn += fabs(n.y);
    } 
#line 301
{
      n.z = val(as.z,0,0,0) - val(as.z,0,0,1);
      nn += fabs(n.z);
    }
    if (nn == 0.)
      val(c,0,0,0) = val(as.x,0,0,0);
    else {
      
 n.x /= nn;
 
#line 309
n.y /= nn;
 
#line 309
n.z /= nn;






      double alpha = 0., ni = 0.;
      for (int i = 0; i <= 1; i++)
 for (int j = 0; j <= 1; j++)
   {
     if (val(p.x,0,i,j) > 0. && val(p.x,0,i,j) < 1.) {
       double a = sign(val(Phi,0,i,j) - val)*(val(p.x,0,i,j) - 0.5);
       alpha += n.x*a + n.y*(i - 0.5) + n.z*(j - 0.5);
       ni++;
     }
     
#line 320
if (val(p.y,j,0,i) > 0. && val(p.y,j,0,i) < 1.) {
       double a = sign(val(Phi,j,0,i) - val)*(val(p.y,j,0,i) - 0.5);
       alpha += n.y*a + n.z*(i - 0.5) + n.x*(j - 0.5);
       ni++;
     }
     
#line 320
if (val(p.z,i,j,0) > 0. && val(p.z,i,j,0) < 1.) {
       double a = sign(val(Phi,i,j,0) - val)*(val(p.z,i,j,0) - 0.5);
       alpha += n.z*a + n.x*(i - 0.5) + n.y*(j - 0.5);
       ni++;
     }}




      if (ni == 0)
 val(c,0,0,0) = val(as.x,0,0,0);
      else if (ni < 3 || ni > 6)
 val(c,0,0,0) = 0.;
      else
 val(c,0,0,0) = plane_volume (n, alpha/ni);
    }
  }end_foreach();}delete((scalar*)((vector[]){p,{{-1},{-1},{-1}}}));
#line 351 "/home/lisergey/basilisk/src/fractions.h"
end_tracing("fractions","/home/lisergey/basilisk/src/fractions.h",351);}
#line 395 "/home/lisergey/basilisk/src/fractions.h"
coord youngs_normal (Point point, scalar c)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;
  coord n;
  double nn = 0.;
  if (!(3 == 2)) qassert ("/home/lisergey/basilisk/src/fractions.h", 399, "dimension == 2");
   {
    n.x = (val(c,-1,1,0) + 2.*val(c,-1,0,0) + val(c,-1,-1,0) -
    val(c,+1,1,0) - 2.*val(c,+1,0,0) - val(c,+1,-1,0));
    nn += fabs(n.x);
  } 
#line 400
{
    n.y = (val(c,0,-1,1) + 2.*val(c,0,-1,0) + val(c,0,-1,-1) -
    val(c,0,+1,1) - 2.*val(c,0,+1,0) - val(c,0,+1,-1));
    nn += fabs(n.y);
  } 
#line 400
{
    n.z = (val(c,1,0,-1) + 2.*val(c,0,0,-1) + val(c,-1,0,-1) -
    val(c,1,0,+1) - 2.*val(c,0,0,+1) - val(c,-1,0,+1));
    nn += fabs(n.z);
  }

  if (nn > 0.)
    {
      n.x /= nn;
      
#line 408
n.y /= nn;
      
#line 408
n.z /= nn;}
  else
    n.x = 1.;
  return n;
}





coord facet_normal (Point point, scalar c, vector s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;
  if (s.x.i >= 0) {
    coord n;
    double nn = 0.;
     {
      n.x = val(s.x,0,0,0) - val(s.x,1,0,0);
      nn += fabs(n.x);
    } 
#line 423
{
      n.y = val(s.y,0,0,0) - val(s.y,0,1,0);
      nn += fabs(n.y);
    } 
#line 423
{
      n.z = val(s.z,0,0,0) - val(s.z,0,0,1);
      nn += fabs(n.z);
    }
    if (nn > 0.)
      {
 n.x /= nn;
 
#line 429
n.y /= nn;
 
#line 429
n.z /= nn;}
    else
      {
 n.x = 1./3;
 
#line 432
n.y = 1./3;
 
#line 432
n.z = 1./3;}
    return n;
  }
  return mycs (point, c);
}






#line 418
static void _stencil_facet_normal (Point point, scalar c, vector s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;
  if (s.x.i >= 0) {    
    
    
     { 
_stencil_val(s.x,0,0,0); _stencil_val(s.x,1,0,0);  
       
       
    
#line 426
} 
#line 423
{ 
_stencil_val(s.y,0,0,0); _stencil_val(s.y,0,1,0);  
       
       
    
#line 426
} 
#line 423
{ 
_stencil_val(s.z,0,0,0); _stencil_val(s.z,0,0,1);  
       
       
    
#line 426
}
      
   
       
    
       
  
    return ;
  } 
_stencil_mycs (point, c);
  
#line 435
return;
}
#line 445 "/home/lisergey/basilisk/src/fractions.h"
     
void reconstruction (const scalar c, vector n, scalar alpha)
{tracing("reconstruction","/home/lisergey/basilisk/src/fractions.h",446);
  foreach_stencil(1,{(NonLocal[]){{"n","vector",(void *)&n,NULL,0},{"alpha","scalar",(void *)&alpha,NULL,0},{"c","scalar",(void *)&c,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"\n\n\n\n\n// #line 62 \"/home/lisergey/basilisk/src/geometry.h\"\nreal plane_alpha (real c, coord n)\n{\n  real alpha;\n  coord n1;\n\n  n1.x = fabs (n.x); n1.y = fabs (n.y); n1.z = fabs (n.z);\n\n  real m1, m2, m3;\n  m1 = min(n1.x, n1.y);\n  m3 = max(n1.x, n1.y);\n  m2 = n1.z;\n  if (m2 < m1) {\n    real tmp = m1;\n    m1 = m2;\n    m2 = tmp;\n  }\n  else if (m2 > m3) {\n    real tmp = m3;\n    m3 = m2;\n    m2 = tmp;\n  }\n  real m12 = m1 + m2;\n  real pr = max(6.*m1*m2*m3, 1e-50);\n  real V1 = m1*m1*m1/pr;\n  real V2 = V1 + (m2 - m1)/(2.*m3), V3;\n  real mm;\n  if (m3 < m12) {\n    mm = m3;\n    V3 = (m3*m3*(3.*m12 - m3) + m1*m1*(m1 - 3.*m3) + m2*m2*(m2 - 3.*m3))/pr;\n  }\n  else {\n    mm = m12;\n    V3 = mm/(2.*m3);\n  }\n\n  c = clamp (c, 0., 1.);\n  real ch = min(c, 1. - c);\n  if (ch < V1)\n    alpha = pow (pr*ch, 1./3.);\n  else if (ch < V2)\n    alpha = (m1 + sqrt(m1*m1 + 8.*m2*m3*(ch - V1)))/2.;\n  else if (ch < V3) {\n    real p12 = sqrt (2.*m1*m2);\n    real q = 3.*(m12 - 2.*m3*ch)/(4.*p12);\n    real teta = acos(clamp(q,-1.,1.))/3.;\n    real cs = cos(teta);\n    alpha = p12*(sqrt(3.*(1. - cs*cs)) - cs) + m12;\n  }\n  else if (m12 <= m3)\n    alpha = m3*ch + mm/2.;\n  else {\n    real p = m1*(m2 + m3) + m2*m3 - 1./4., p12 = sqrt(p);\n    real q = 3.*m1*m2*m3*(1./2. - ch)/(2.*p*p12);\n    real teta = acos(clamp(q,-1.,1.))/3.;\n    real cs = cos(teta);\n    alpha = p12*(sqrt(3.*(1. - cs*cs)) - cs) + 1./2.;\n  }\n  if (c > 1./2.) alpha = 1. - alpha;\n\n  if (n.x < 0.)\n    alpha += n.x;\n  if (n.y < 0.)\n    alpha += n.y;\n  if (n.z < 0.)\n    alpha += n.z;\n\n  return alpha - (n.x + n.y + n.z)/2.;;\n}\n// # 13 \"/home/lisergey/basilisk/src/fractions.h\" 2\n\n\n\n\n\n\n\n// # 1 \"myc.h\" 1\n// # 1 \"/home/lisergey/basilisk/src/myc.h\"\n// # 16 \"/home/lisergey/basilisk/src/myc.h\"\ncoord mycs (Point point, scalar c)\n{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;\n  real m1,m2,m[4][3],t0,t1,t2;\n  int cn;\n\n\n\n  m1 = val(c,-1,0,-1) + val(c,-1,0,1) + val(c,-1,-1,0) + val(c,-1,1,0) +\n       val(c,-1,0,0);\n  m2 = val(c,1,0,-1) + val(c,1,0,1) + val(c,1,-1,0) + val(c,1,1,0) +\n       val(c,1,0,0);\n  m[0][0] = m1 > m2 ? 1. : -1.;\n\n  m1 = val(c,-1,-1,0)+ val(c,1,-1,0)+ val(c,0,-1,0);\n  m2 = val(c,-1,1,0)+ val(c,1,1,0)+ val(c,0,1,0);\n  m[0][1] = 0.5*(m1-m2);\n\n  m1 = val(c,-1,0,-1)+ val(c,1,0,-1)+ val(c,0,0,-1);\n  m2 = val(c,-1,0,1)+ val(c,1,0,1)+ val(c,0,0,1);\n  m[0][2] = 0.5*(m1-m2);\n\n\n\n  m1 = val(c,-1,-1,0) + val(c,-1,1,0) + val(c,-1,0,0);\n  m2 = val(c,1,-1,0) + val(c,1,1,0) + val(c,1,0,0);\n  m[1][0] = 0.5*(m1-m2);\n\n  m1 = val(c,0,-1,-1) + val(c,0,-1,1) + val(c,1,-1,0) + val(c,-1,-1,0) +\n       val(c,0,-1,0);\n  m2 = val(c,0,1,-1) + val(c,0,1,1) + val(c,1,1,0) + val(c,-1,1,0) +\n       val(c,0,1,0);\n  m[1][1] = m1 > m2 ? 1. : -1.;\n\n  m1 = val(c,0,-1,-1)+ val(c,0,0,-1)+ val(c,0,1,-1);\n  m2 = val(c,0,-1,1)+ val(c,0,0,1)+ val(c,0,1,1);\n  m[1][2] = 0.5*(m1-m2);\n\n\n\n\n  m1 = val(c,-1,0,-1)+ val(c,-1,0,1)+ val(c,-1,0,0);\n  m2 = val(c,1,0,-1)+ val(c,1,0,1)+ val(c,1,0,0);\n  m[2][0] = 0.5*(m1-m2);\n\n  m1 = val(c,0,-1,-1)+ val(c,0,-1,1)+ val(c,0,-1,0);\n  m2 = val(c,0,1,-1)+ val(c,0,1,1)+ val(c,0,1,0);\n  m[2][1] = 0.5*(m1-m2);\n\n  m1 = val(c,-1,0,-1) + val(c,1,0,-1) + val(c,0,-1,-1) + val(c,0,1,-1) +\n       val(c,0,0,-1);\n  m2 = val(c,-1,0,1) + val(c,1,0,1) + val(c,0,-1,1) + val(c,0,1,1) +\n       val(c,0,0,1);\n  m[2][2] = m1 > m2 ? 1. : -1.;\n\n\n  t0 = fabs(m[0][0]) + fabs(m[0][1]) + fabs(m[0][2]);\n  m[0][0] /= t0;\n  m[0][1] /= t0;\n  m[0][2] /= t0;\n\n  t0 = fabs(m[1][0]) + fabs(m[1][1]) + fabs(m[1][2]);\n  m[1][0] /= t0;\n  m[1][1] /= t0;\n  m[1][2] /= t0;\n\n  t0 = fabs(m[2][0]) + fabs(m[2][1]) + fabs(m[2][2]);\n  m[2][0] /= t0;\n  m[2][1] /= t0;\n  m[2][2] /= t0;\n\n\n  t0 = fabs(m[0][0]);\n  t1 = fabs(m[1][1]);\n  t2 = fabs(m[2][2]);\n\n  cn = 0;\n  if (t1 > t0) {\n    t0 = t1;\n    cn = 1;\n  }\n  if (t2 > t0)\n    cn = 2;\n\n\n  m1 = val(c,-1,-1,-1) + val(c,-1,1,-1) + val(c,-1,-1,1) + val(c,-1,1,1) +\n       2.*(val(c,-1,-1,0) + val(c,-1,1,0) + val(c,-1,0,-1) + val(c,-1,0,1)) +\n       4.*val(c,-1,0,0);\n  m2 = val(c,1,-1,-1) + val(c,1,1,-1) + val(c,1,-1,1) + val(c,1,1,1) +\n       2.*(val(c,1,-1,0) + val(c,1,1,0) + val(c,1,0,-1) + val(c,1,0,1)) +\n       4.*val(c,1,0,0);\n  m[3][0] = m1 - m2;\n\n  m1 = val(c,-1,-1,-1) + val(c,-1,-1,1) + val(c,1,-1,-1) + val(c,1,-1,1) +\n       2.*( val(c,-1,-1,0) + val(c,1,-1,0) + val(c,0,-1,-1) + val(c,0,-1,1)) +\n       4.*val(c,0,-1,0);\n  m2 = val(c,-1,1,-1) + val(c,-1,1,1) + val(c,1,1,-1) + val(c,1,1,1) +\n       2.*(val(c,-1,1,0) + val(c,1,1,0) + val(c,0,1,-1) + val(c,0,1,1)) +\n       4.*val(c,0,1,0);\n  m[3][1] = m1 - m2;\n\n  m1 = val(c,-1,-1,-1) + val(c,-1,1,-1) + val(c,1,-1,-1) + val(c,1,1,-1) +\n       2.*(val(c,-1,0,-1) + val(c,1,0,-1) + val(c,0,-1,-1) + val(c,0,1,-1)) +\n       4.*val(c,0,0,-1);\n  m2 = val(c,-1,-1,1) + val(c,-1,1,1) + val(c,1,-1,1) + val(c,1,1,1) +\n       2.*(val(c,-1,0,1) + val(c,1,0,1) + val(c,0,-1,1) + val(c,0,1,1)) +\n       4.*val(c,0,0,1);\n  m[3][2] = m1 - m2;\n\n\n  t0 = fabs(m[3][0]) + fabs(m[3][1]) + fabs(m[3][2]);\n  if (t0 < 1e-30) {\n    coord mxyz = {1., 0., 0.};\n    return mxyz;\n  }\n\n  m[3][0] /= t0;\n  m[3][1] /= t0;\n  m[3][2] /= t0;\n\n\n  t0 = fabs (m[3][0]);\n  t1 = fabs (m[3][1]);\n  t2 = fabs (m[3][2]);\n  if (t1 > t0)\n    t0 = t1;\n  if (t2 > t0)\n    t0 = t2;\n\n  if (fabs(m[cn][cn]) > t0)\n    cn = 3;\n\n\n  coord mxyz = {m[cn][0], m[cn][1], m[cn][2]};\n  return mxyz;\n}"," \n// #line 448 \"/home/lisergey/basilisk/src/fractions.h\"\n{\n\n\n\n\n\n    if (val(c,0,0,0) <= 0. || val(c,0,0,0) >= 1.) {\n      val_out_(alpha,0,0,0) = 0.;\n      \n val_out_(n.x,0,0,0) = 0.;\n \n// #line 457\nval_out_(n.y,0,0,0) = 0.;\n \n// #line 457\nval_out_(n.z,0,0,0) = 0.;\n    }\n    else {\n\n\n\n\n\n\n      coord m = mycs (point, c);\n      \n val_out_(n.x,0,0,0) = m.x;\n \n// #line 468\nval_out_(n.y,0,0,0) = m.y;\n \n// #line 468\nval_out_(n.z,0,0,0) = m.z;\n      val_out_(alpha,0,0,0) = plane_alpha (val(c,0,0,0), m);\n    }\n  }"}) {





_stencil_val(c,0,0,0); _stencil_val(c,0,0,0);{ {
      _stencil_val_a(alpha,0,0,0);  
      
 {_stencil_val_a(n.x,0,0,0);  }
 
#line 457
{_stencil_val_a(n.y,0,0,0);  }
 
#line 457
{_stencil_val_a(n.z,0,0,0);  }
    } 
{  






       _stencil_mycs (point, c);
      
 {_stencil_val_a(n.x,0,0,0);  }
 
#line 468
{_stencil_val_a(n.y,0,0,0);  }
 
#line 468
{_stencil_val_a(n.z,0,0,0);  }
_stencil_val(c,0,0,0);
      
#line 469
_stencil_val_a(alpha,0,0,0);    
    }}





          
    
  
#line 471
}end_foreach_stencil();
  {
#line 448
foreach() {





    if (val(c,0,0,0) <= 0. || val(c,0,0,0) >= 1.) {
      val(alpha,0,0,0) = 0.;
      
 val(n.x,0,0,0) = 0.;
 
#line 457
val(n.y,0,0,0) = 0.;
 
#line 457
val(n.z,0,0,0) = 0.;
    }
    else {






      coord m = mycs (point, c);
      
 val(n.x,0,0,0) = m.x;
 
#line 468
val(n.y,0,0,0) = m.y;
 
#line 468
val(n.z,0,0,0) = m.z;
      val(alpha,0,0,0) = plane_alpha (val(c,0,0,0), m);
    }
  }end_foreach();}
#line 480 "/home/lisergey/basilisk/src/fractions.h"
  
    _attribute[n.x.i].refine = _attribute[n.x.i].prolongation = refine_injection;
    
#line 481
_attribute[n.y.i].refine = _attribute[n.y.i].prolongation = refine_injection;
    
#line 481
_attribute[n.z.i].refine = _attribute[n.z.i].prolongation = refine_injection;




  _attribute[alpha.i].n = n;
  _attribute[alpha.i].refine = _attribute[alpha.i].prolongation = alpha_refine;

end_tracing("reconstruction","/home/lisergey/basilisk/src/fractions.h",489);}
#line 509 "/home/lisergey/basilisk/src/fractions.h"
     
void output_facets (scalar c, FILE * fp, vector s)
{tracing("output_facets","/home/lisergey/basilisk/src/fractions.h",510);
  foreach_stencil(1,{(NonLocal[]){{"fp","not implemented yet",(void *)fp,NULL,1},{"cube_connect","int",(void *)cube_connect,(int[]){1224,0},0},{"cube_edge","coord",(void *)cube_edge,(int[]){122,0},0},{"s","vector",(void *)&s,NULL,0},{"c","scalar",(void *)&c,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"\n\n\n// #line 334 \"/home/lisergey/basilisk/src/geometry.h\"\nint facets (coord n, real alpha, coord v[12], real h)\n{\n  coord a[12];\n  int orient[12];\n\n  for (int i = 0; i < 12; i++) {\n    coord e, d;\n    real den = 0., t = alpha;\n     {\n      d.x = h*(cube_edge[i][0].x - 0.5);\n      e.x = h*(cube_edge[i][1].x - 0.5);\n      den += n.x*(e.x - d.x);\n      t -= n.x*d.x;\n    } \n// #line 342\n{\n      d.y = h*(cube_edge[i][0].y - 0.5);\n      e.y = h*(cube_edge[i][1].y - 0.5);\n      den += n.y*(e.y - d.y);\n      t -= n.y*d.y;\n    } \n// #line 342\n{\n      d.z = h*(cube_edge[i][0].z - 0.5);\n      e.z = h*(cube_edge[i][1].z - 0.5);\n      den += n.z*(e.z - d.z);\n      t -= n.z*d.z;\n    }\n    orient[i] = -1;\n    if (fabs (den) > 1e-10) {\n      t /= den;\n      if (t >= 0. && t < 1.) {\n real s = - alpha;\n  {\n   a[i].x = d.x + t*(e.x - d.x);\n   s += n.x*e.x;\n } \n// #line 353\n{\n   a[i].y = d.y + t*(e.y - d.y);\n   s += n.y*e.y;\n } \n// #line 353\n{\n   a[i].z = d.z + t*(e.z - d.z);\n   s += n.z*e.z;\n }\n orient[i] = (s > 0.);\n      }\n    }\n  }\n\n  for (int i = 0; i < 12; i++) {\n    int nv = 0, e = i;\n    while (orient[e] >= 0) {\n      int m = 0, * ne = cube_connect[e][orient[e]];\n      v[nv++] = a[e];\n      orient[e] = -1;\n      while (m < 3 && orient[e] < 0)\n e = ne[m++];\n    }\n    if (nv > 2)\n      return nv;\n  }\n  return 0;\n}\n\n\n\n\n// #line 62 \"/home/lisergey/basilisk/src/geometry.h\"\nreal plane_alpha (real c, coord n)\n{\n  real alpha;\n  coord n1;\n\n  n1.x = fabs (n.x); n1.y = fabs (n.y); n1.z = fabs (n.z);\n\n  real m1, m2, m3;\n  m1 = min(n1.x, n1.y);\n  m3 = max(n1.x, n1.y);\n  m2 = n1.z;\n  if (m2 < m1) {\n    real tmp = m1;\n    m1 = m2;\n    m2 = tmp;\n  }\n  else if (m2 > m3) {\n    real tmp = m3;\n    m3 = m2;\n    m2 = tmp;\n  }\n  real m12 = m1 + m2;\n  real pr = max(6.*m1*m2*m3, 1e-50);\n  real V1 = m1*m1*m1/pr;\n  real V2 = V1 + (m2 - m1)/(2.*m3), V3;\n  real mm;\n  if (m3 < m12) {\n    mm = m3;\n    V3 = (m3*m3*(3.*m12 - m3) + m1*m1*(m1 - 3.*m3) + m2*m2*(m2 - 3.*m3))/pr;\n  }\n  else {\n    mm = m12;\n    V3 = mm/(2.*m3);\n  }\n\n  c = clamp (c, 0., 1.);\n  real ch = min(c, 1. - c);\n  if (ch < V1)\n    alpha = pow (pr*ch, 1./3.);\n  else if (ch < V2)\n    alpha = (m1 + sqrt(m1*m1 + 8.*m2*m3*(ch - V1)))/2.;\n  else if (ch < V3) {\n    real p12 = sqrt (2.*m1*m2);\n    real q = 3.*(m12 - 2.*m3*ch)/(4.*p12);\n    real teta = acos(clamp(q,-1.,1.))/3.;\n    real cs = cos(teta);\n    alpha = p12*(sqrt(3.*(1. - cs*cs)) - cs) + m12;\n  }\n  else if (m12 <= m3)\n    alpha = m3*ch + mm/2.;\n  else {\n    real p = m1*(m2 + m3) + m2*m3 - 1./4., p12 = sqrt(p);\n    real q = 3.*m1*m2*m3*(1./2. - ch)/(2.*p*p12);\n    real teta = acos(clamp(q,-1.,1.))/3.;\n    real cs = cos(teta);\n    alpha = p12*(sqrt(3.*(1. - cs*cs)) - cs) + 1./2.;\n  }\n  if (c > 1./2.) alpha = 1. - alpha;\n\n  if (n.x < 0.)\n    alpha += n.x;\n  if (n.y < 0.)\n    alpha += n.y;\n  if (n.z < 0.)\n    alpha += n.z;\n\n  return alpha - (n.x + n.y + n.z)/2.;;\n}\n// # 13 \"/home/lisergey/basilisk/src/fractions.h\" 2\n\n\n\n\n\n\n\n// # 1 \"myc.h\" 1\n// # 1 \"/home/lisergey/basilisk/src/myc.h\"\n// # 16 \"/home/lisergey/basilisk/src/myc.h\"\ncoord mycs (Point point, scalar c)\n{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;\n  real m1,m2,m[4][3],t0,t1,t2;\n  int cn;\n\n\n\n  m1 = val(c,-1,0,-1) + val(c,-1,0,1) + val(c,-1,-1,0) + val(c,-1,1,0) +\n       val(c,-1,0,0);\n  m2 = val(c,1,0,-1) + val(c,1,0,1) + val(c,1,-1,0) + val(c,1,1,0) +\n       val(c,1,0,0);\n  m[0][0] = m1 > m2 ? 1. : -1.;\n\n  m1 = val(c,-1,-1,0)+ val(c,1,-1,0)+ val(c,0,-1,0);\n  m2 = val(c,-1,1,0)+ val(c,1,1,0)+ val(c,0,1,0);\n  m[0][1] = 0.5*(m1-m2);\n\n  m1 = val(c,-1,0,-1)+ val(c,1,0,-1)+ val(c,0,0,-1);\n  m2 = val(c,-1,0,1)+ val(c,1,0,1)+ val(c,0,0,1);\n  m[0][2] = 0.5*(m1-m2);\n\n\n\n  m1 = val(c,-1,-1,0) + val(c,-1,1,0) + val(c,-1,0,0);\n  m2 = val(c,1,-1,0) + val(c,1,1,0) + val(c,1,0,0);\n  m[1][0] = 0.5*(m1-m2);\n\n  m1 = val(c,0,-1,-1) + val(c,0,-1,1) + val(c,1,-1,0) + val(c,-1,-1,0) +\n       val(c,0,-1,0);\n  m2 = val(c,0,1,-1) + val(c,0,1,1) + val(c,1,1,0) + val(c,-1,1,0) +\n       val(c,0,1,0);\n  m[1][1] = m1 > m2 ? 1. : -1.;\n\n  m1 = val(c,0,-1,-1)+ val(c,0,0,-1)+ val(c,0,1,-1);\n  m2 = val(c,0,-1,1)+ val(c,0,0,1)+ val(c,0,1,1);\n  m[1][2] = 0.5*(m1-m2);\n\n\n\n\n  m1 = val(c,-1,0,-1)+ val(c,-1,0,1)+ val(c,-1,0,0);\n  m2 = val(c,1,0,-1)+ val(c,1,0,1)+ val(c,1,0,0);\n  m[2][0] = 0.5*(m1-m2);\n\n  m1 = val(c,0,-1,-1)+ val(c,0,-1,1)+ val(c,0,-1,0);\n  m2 = val(c,0,1,-1)+ val(c,0,1,1)+ val(c,0,1,0);\n  m[2][1] = 0.5*(m1-m2);\n\n  m1 = val(c,-1,0,-1) + val(c,1,0,-1) + val(c,0,-1,-1) + val(c,0,1,-1) +\n       val(c,0,0,-1);\n  m2 = val(c,-1,0,1) + val(c,1,0,1) + val(c,0,-1,1) + val(c,0,1,1) +\n       val(c,0,0,1);\n  m[2][2] = m1 > m2 ? 1. : -1.;\n\n\n  t0 = fabs(m[0][0]) + fabs(m[0][1]) + fabs(m[0][2]);\n  m[0][0] /= t0;\n  m[0][1] /= t0;\n  m[0][2] /= t0;\n\n  t0 = fabs(m[1][0]) + fabs(m[1][1]) + fabs(m[1][2]);\n  m[1][0] /= t0;\n  m[1][1] /= t0;\n  m[1][2] /= t0;\n\n  t0 = fabs(m[2][0]) + fabs(m[2][1]) + fabs(m[2][2]);\n  m[2][0] /= t0;\n  m[2][1] /= t0;\n  m[2][2] /= t0;\n\n\n  t0 = fabs(m[0][0]);\n  t1 = fabs(m[1][1]);\n  t2 = fabs(m[2][2]);\n\n  cn = 0;\n  if (t1 > t0) {\n    t0 = t1;\n    cn = 1;\n  }\n  if (t2 > t0)\n    cn = 2;\n\n\n  m1 = val(c,-1,-1,-1) + val(c,-1,1,-1) + val(c,-1,-1,1) + val(c,-1,1,1) +\n       2.*(val(c,-1,-1,0) + val(c,-1,1,0) + val(c,-1,0,-1) + val(c,-1,0,1)) +\n       4.*val(c,-1,0,0);\n  m2 = val(c,1,-1,-1) + val(c,1,1,-1) + val(c,1,-1,1) + val(c,1,1,1) +\n       2.*(val(c,1,-1,0) + val(c,1,1,0) + val(c,1,0,-1) + val(c,1,0,1)) +\n       4.*val(c,1,0,0);\n  m[3][0] = m1 - m2;\n\n  m1 = val(c,-1,-1,-1) + val(c,-1,-1,1) + val(c,1,-1,-1) + val(c,1,-1,1) +\n       2.*( val(c,-1,-1,0) + val(c,1,-1,0) + val(c,0,-1,-1) + val(c,0,-1,1)) +\n       4.*val(c,0,-1,0);\n  m2 = val(c,-1,1,-1) + val(c,-1,1,1) + val(c,1,1,-1) + val(c,1,1,1) +\n       2.*(val(c,-1,1,0) + val(c,1,1,0) + val(c,0,1,-1) + val(c,0,1,1)) +\n       4.*val(c,0,1,0);\n  m[3][1] = m1 - m2;\n\n  m1 = val(c,-1,-1,-1) + val(c,-1,1,-1) + val(c,1,-1,-1) + val(c,1,1,-1) +\n       2.*(val(c,-1,0,-1) + val(c,1,0,-1) + val(c,0,-1,-1) + val(c,0,1,-1)) +\n       4.*val(c,0,0,-1);\n  m2 = val(c,-1,-1,1) + val(c,-1,1,1) + val(c,1,-1,1) + val(c,1,1,1) +\n       2.*(val(c,-1,0,1) + val(c,1,0,1) + val(c,0,-1,1) + val(c,0,1,1)) +\n       4.*val(c,0,0,1);\n  m[3][2] = m1 - m2;\n\n\n  t0 = fabs(m[3][0]) + fabs(m[3][1]) + fabs(m[3][2]);\n  if (t0 < 1e-30) {\n    coord mxyz = {1., 0., 0.};\n    return mxyz;\n  }\n\n  m[3][0] /= t0;\n  m[3][1] /= t0;\n  m[3][2] /= t0;\n\n\n  t0 = fabs (m[3][0]);\n  t1 = fabs (m[3][1]);\n  t2 = fabs (m[3][2]);\n  if (t1 > t0)\n    t0 = t1;\n  if (t2 > t0)\n    t0 = t2;\n\n  if (fabs(m[cn][cn]) > t0)\n    cn = 3;\n\n\n  coord mxyz = {m[cn][0], m[cn][1], m[cn][2]};\n  return mxyz;\n}\n\n\n\n\n\n\n// #line 418 \"/home/lisergey/basilisk/src/fractions.h\"\ncoord facet_normal (Point point, scalar c, vector s)\n{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;\n  if (s.x.i >= 0) {\n    coord n;\n    real nn = 0.;\n     {\n      n.x = val(s.x,0,0,0) - val(s.x,1,0,0);\n      nn += fabs(n.x);\n    } \n// #line 423\n{\n      n.y = val(s.y,0,0,0) - val(s.y,0,1,0);\n      nn += fabs(n.y);\n    } \n// #line 423\n{\n      n.z = val(s.z,0,0,0) - val(s.z,0,0,1);\n      nn += fabs(n.z);\n    }\n    if (nn > 0.)\n      {\n n.x /= nn;\n \n// #line 429\nn.y /= nn;\n \n// #line 429\nn.z /= nn;}\n    else\n      {\n n.x = 1./3;\n \n// #line 432\nn.y = 1./3;\n \n// #line 432\nn.z = 1./3;}\n    return n;\n  }\n  return mycs (point, c);\n}","\n    \n// #line 513 \"/home/lisergey/basilisk/src/fractions.h\"\nif (val(c,0,0,0) > 1e-6 && val(c,0,0,0) < 1. - 1e-6) {\n      coord n = facet_normal (point, c, s);\n      real alpha = plane_alpha (val(c,0,0,0), n);\n// # 525 \"/home/lisergey/basilisk/src/fractions.h\"\n      coord v[12];\n      int m = facets (n, alpha, v, 1.);\n      for (int i = 0; i < m; i++)\n fprintf (fp, \"%g %g %g\\n\",\n   x + v[i].x*Delta, y + v[i].y*Delta, z + v[i].z*Delta);\n      if (m > 0)\n fputc ('\\n', fp);\n\n    }"})
    {_stencil_val(c,0,0,0); _stencil_val(c,0,0,0); {  
       _stencil_facet_normal (point, c, s);     
      _stencil_val(c,0,0,0);        
#line 525 "/home/lisergey/basilisk/src/fractions.h"
      
                 
              
     
 
        
 

    }        }end_foreach_stencil();
  {
#line 512
foreach()
    if (val(c,0,0,0) > 1e-6 && val(c,0,0,0) < 1. - 1e-6) {
      coord n = facet_normal (point, c, s);
      double alpha = plane_alpha (val(c,0,0,0), n);
#line 525 "/home/lisergey/basilisk/src/fractions.h"
      coord v[12];
      int m = facets (n, alpha, v, 1.);
      for (int i = 0; i < m; i++)
 fprintf (fp, "%g %g %g\n",
   x + v[i].x*Delta, y + v[i].y*Delta, z + v[i].z*Delta);
      if (m > 0)
 fputc ('\n', fp);

    }end_foreach();}

  fflush (fp);
end_tracing("output_facets","/home/lisergey/basilisk/src/fractions.h",536);}







     
double interface_area (scalar c)
{tracing("interface_area","/home/lisergey/basilisk/src/fractions.h",545);
  double area = 0.;
  foreach_stencil (1,{(NonLocal[]){{"area","double",(void *)&area,NULL,0,'+'},{"c","scalar",(void *)&c,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"\n\n\n\n\n\n\n\n// #line 382 \"/home/lisergey/basilisk/src/geometry.h\"\nreal line_length_center (coord m, real alpha, coord * p)\n{\n  alpha += (m.x + m.y)/2.;\n\n  coord n = m;\n  \n    if (n.x < 0.) {\n      alpha -= n.x;\n      n.x = - n.x;\n    }\n    \n// #line 388\nif (n.y < 0.) {\n      alpha -= n.y;\n      n.y = - n.y;\n    }\n\n  p->x = p->y = p->z = 0.;\n\n  if (alpha <= 0. || alpha >= n.x + n.y)\n    return 0.;\n\n  \n    if (n.x < 1e-4) {\n      p->x = 0.;\n      p->y = (m.y < 0. ? 1. - alpha : alpha) - 0.5;\n      return 1.;\n    }\n    \n// #line 399\nif (n.y < 1e-4) {\n      p->y = 0.;\n      p->x = (m.x < 0. ? 1. - alpha : alpha) - 0.5;\n      return 1.;\n    }\n\n  if (alpha >= n.x) {\n    p->x += 1.;\n    p->y += (alpha - n.x)/n.y;\n  }\n  else\n    p->x += alpha/n.x;\n\n  real ax = p->x, ay = p->y;\n  if (alpha >= n.y) {\n    p->y += 1.;\n    ay -= 1.;\n    p->x += (alpha - n.y)/n.x;\n    ax -= (alpha - n.y)/n.x;\n  }\n  else {\n    p->y += alpha/n.y;\n    ay -= alpha/n.y;\n  }\n\n   {\n    p->x /= 2.;\n    p->x = clamp (p->x, 0., 1.);\n    if (m.x < 0.)\n      p->x = 1. - p->x;\n    p->x -= 0.5;\n  } \n// #line 424\n{\n    p->y /= 2.;\n    p->y = clamp (p->y, 0., 1.);\n    if (m.y < 0.)\n      p->y = 1. - p->y;\n    p->y -= 0.5;\n  }\n\n  return sqrt (ax*ax + ay*ay);\n}\n\n\n\n\n\n// #line 438 \"/home/lisergey/basilisk/src/geometry.h\"\nreal plane_area_center (coord m, real alpha, coord * p)\n{\n  \n    if (fabs (m.x) < 1e-4) {\n      coord n, q;\n      ((real *)&n)[0] = m.y;\n      ((real *)&n)[1] = m.z;\n      real length = line_length_center (n, alpha, &q);\n      p->x = 0.;\n      p->y = ((real *)&q)[0];\n      p->z = ((real *)&q)[1];\n      return length;\n    }\n    \n// #line 441\nif (fabs (m.y) < 1e-4) {\n      coord n, q;\n      ((real *)&n)[0] = m.z;\n      ((real *)&n)[1] = m.x;\n      real length = line_length_center (n, alpha, &q);\n      p->y = 0.;\n      p->z = ((real *)&q)[0];\n      p->x = ((real *)&q)[1];\n      return length;\n    }\n    \n// #line 441\nif (fabs (m.z) < 1e-4) {\n      coord n, q;\n      ((real *)&n)[0] = m.x;\n      ((real *)&n)[1] = m.y;\n      real length = line_length_center (n, alpha, &q);\n      p->z = 0.;\n      p->x = ((real *)&q)[0];\n      p->y = ((real *)&q)[1];\n      return length;\n    }\n\n  alpha += (m.x + m.y + m.z)/2.;\n  coord n = m;\n  \n    if (n.x < 0.) {\n      alpha -= n.x;\n      n.x = - n.x;\n    }\n    \n// #line 455\nif (n.y < 0.) {\n      alpha -= n.y;\n      n.y = - n.y;\n    }\n    \n// #line 455\nif (n.z < 0.) {\n      alpha -= n.z;\n      n.z = - n.z;\n    }\n\n  real amax = n.x + n.y + n.z;\n  if (alpha < 0. || alpha > amax) {\n    p->x = p->y = p->z = 0.;\n    return 0.;\n  }\n\n  real area = sq(alpha);\n  p->x = p->y = p->z = area*alpha;\n\n   {\n    real b = alpha - n.x;\n    if (b > 0.) {\n      area -= b*b;\n      p->x -= b*b*(2.*n.x + alpha);\n      p->y -= b*b*b;\n      p->z -= b*b*b;\n    }\n  } \n// #line 469\n{\n    real b = alpha - n.y;\n    if (b > 0.) {\n      area -= b*b;\n      p->y -= b*b*(2.*n.y + alpha);\n      p->z -= b*b*b;\n      p->x -= b*b*b;\n    }\n  } \n// #line 469\n{\n    real b = alpha - n.z;\n    if (b > 0.) {\n      area -= b*b;\n      p->z -= b*b*(2.*n.z + alpha);\n      p->x -= b*b*b;\n      p->y -= b*b*b;\n    }\n  }\n\n  amax = alpha - amax;\n   {\n    real b = amax + n.x;\n    if (b > 0.) {\n      area += b*b;\n      p->y += b*b*(2.*n.y + alpha - n.z);\n      p->z += b*b*(2.*n.z + alpha - n.y);\n      p->x += b*b*b;\n    }\n  } \n// #line 480\n{\n    real b = amax + n.y;\n    if (b > 0.) {\n      area += b*b;\n      p->z += b*b*(2.*n.z + alpha - n.x);\n      p->x += b*b*(2.*n.x + alpha - n.z);\n      p->y += b*b*b;\n    }\n  } \n// #line 480\n{\n    real b = amax + n.z;\n    if (b > 0.) {\n      area += b*b;\n      p->x += b*b*(2.*n.x + alpha - n.y);\n      p->y += b*b*(2.*n.y + alpha - n.x);\n      p->z += b*b*b;\n    }\n  }\n\n  area *= 3.;\n   {\n    if (area) {\n      p->x /= area*n.x;\n      p->x = clamp (p->x, 0., 1.);\n    }\n    else\n      p->x = 0.;\n    if (m.x < 0.) p->x = 1. - p->x;\n    p->x -= 0.5;\n  } \n// #line 491\n{\n    if (area) {\n      p->y /= area*n.y;\n      p->y = clamp (p->y, 0., 1.);\n    }\n    else\n      p->y = 0.;\n    if (m.y < 0.) p->y = 1. - p->y;\n    p->y -= 0.5;\n  } \n// #line 491\n{\n    if (area) {\n      p->z /= area*n.z;\n      p->z = clamp (p->z, 0., 1.);\n    }\n    else\n      p->z = 0.;\n    if (m.z < 0.) p->z = 1. - p->z;\n    p->z -= 0.5;\n  }\n\n  return area*sqrt (1./(sq(n.x)*sq(n.y)) +\n      1./(sq(n.x)*sq(n.z)) +\n      1./(sq(n.z)*sq(n.y)))/6.;\n}\n\n\n\n\n// #line 62 \"/home/lisergey/basilisk/src/geometry.h\"\nreal plane_alpha (real c, coord n)\n{\n  real alpha;\n  coord n1;\n\n  n1.x = fabs (n.x); n1.y = fabs (n.y); n1.z = fabs (n.z);\n\n  real m1, m2, m3;\n  m1 = min(n1.x, n1.y);\n  m3 = max(n1.x, n1.y);\n  m2 = n1.z;\n  if (m2 < m1) {\n    real tmp = m1;\n    m1 = m2;\n    m2 = tmp;\n  }\n  else if (m2 > m3) {\n    real tmp = m3;\n    m3 = m2;\n    m2 = tmp;\n  }\n  real m12 = m1 + m2;\n  real pr = max(6.*m1*m2*m3, 1e-50);\n  real V1 = m1*m1*m1/pr;\n  real V2 = V1 + (m2 - m1)/(2.*m3), V3;\n  real mm;\n  if (m3 < m12) {\n    mm = m3;\n    V3 = (m3*m3*(3.*m12 - m3) + m1*m1*(m1 - 3.*m3) + m2*m2*(m2 - 3.*m3))/pr;\n  }\n  else {\n    mm = m12;\n    V3 = mm/(2.*m3);\n  }\n\n  c = clamp (c, 0., 1.);\n  real ch = min(c, 1. - c);\n  if (ch < V1)\n    alpha = pow (pr*ch, 1./3.);\n  else if (ch < V2)\n    alpha = (m1 + sqrt(m1*m1 + 8.*m2*m3*(ch - V1)))/2.;\n  else if (ch < V3) {\n    real p12 = sqrt (2.*m1*m2);\n    real q = 3.*(m12 - 2.*m3*ch)/(4.*p12);\n    real teta = acos(clamp(q,-1.,1.))/3.;\n    real cs = cos(teta);\n    alpha = p12*(sqrt(3.*(1. - cs*cs)) - cs) + m12;\n  }\n  else if (m12 <= m3)\n    alpha = m3*ch + mm/2.;\n  else {\n    real p = m1*(m2 + m3) + m2*m3 - 1./4., p12 = sqrt(p);\n    real q = 3.*m1*m2*m3*(1./2. - ch)/(2.*p*p12);\n    real teta = acos(clamp(q,-1.,1.))/3.;\n    real cs = cos(teta);\n    alpha = p12*(sqrt(3.*(1. - cs*cs)) - cs) + 1./2.;\n  }\n  if (c > 1./2.) alpha = 1. - alpha;\n\n  if (n.x < 0.)\n    alpha += n.x;\n  if (n.y < 0.)\n    alpha += n.y;\n  if (n.z < 0.)\n    alpha += n.z;\n\n  return alpha - (n.x + n.y + n.z)/2.;;\n}\n// # 13 \"/home/lisergey/basilisk/src/fractions.h\" 2\n\n\n\n\n\n\n\n// # 1 \"myc.h\" 1\n// # 1 \"/home/lisergey/basilisk/src/myc.h\"\n// # 16 \"/home/lisergey/basilisk/src/myc.h\"\ncoord mycs (Point point, scalar c)\n{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;\n  real m1,m2,m[4][3],t0,t1,t2;\n  int cn;\n\n\n\n  m1 = val(c,-1,0,-1) + val(c,-1,0,1) + val(c,-1,-1,0) + val(c,-1,1,0) +\n       val(c,-1,0,0);\n  m2 = val(c,1,0,-1) + val(c,1,0,1) + val(c,1,-1,0) + val(c,1,1,0) +\n       val(c,1,0,0);\n  m[0][0] = m1 > m2 ? 1. : -1.;\n\n  m1 = val(c,-1,-1,0)+ val(c,1,-1,0)+ val(c,0,-1,0);\n  m2 = val(c,-1,1,0)+ val(c,1,1,0)+ val(c,0,1,0);\n  m[0][1] = 0.5*(m1-m2);\n\n  m1 = val(c,-1,0,-1)+ val(c,1,0,-1)+ val(c,0,0,-1);\n  m2 = val(c,-1,0,1)+ val(c,1,0,1)+ val(c,0,0,1);\n  m[0][2] = 0.5*(m1-m2);\n\n\n\n  m1 = val(c,-1,-1,0) + val(c,-1,1,0) + val(c,-1,0,0);\n  m2 = val(c,1,-1,0) + val(c,1,1,0) + val(c,1,0,0);\n  m[1][0] = 0.5*(m1-m2);\n\n  m1 = val(c,0,-1,-1) + val(c,0,-1,1) + val(c,1,-1,0) + val(c,-1,-1,0) +\n       val(c,0,-1,0);\n  m2 = val(c,0,1,-1) + val(c,0,1,1) + val(c,1,1,0) + val(c,-1,1,0) +\n       val(c,0,1,0);\n  m[1][1] = m1 > m2 ? 1. : -1.;\n\n  m1 = val(c,0,-1,-1)+ val(c,0,0,-1)+ val(c,0,1,-1);\n  m2 = val(c,0,-1,1)+ val(c,0,0,1)+ val(c,0,1,1);\n  m[1][2] = 0.5*(m1-m2);\n\n\n\n\n  m1 = val(c,-1,0,-1)+ val(c,-1,0,1)+ val(c,-1,0,0);\n  m2 = val(c,1,0,-1)+ val(c,1,0,1)+ val(c,1,0,0);\n  m[2][0] = 0.5*(m1-m2);\n\n  m1 = val(c,0,-1,-1)+ val(c,0,-1,1)+ val(c,0,-1,0);\n  m2 = val(c,0,1,-1)+ val(c,0,1,1)+ val(c,0,1,0);\n  m[2][1] = 0.5*(m1-m2);\n\n  m1 = val(c,-1,0,-1) + val(c,1,0,-1) + val(c,0,-1,-1) + val(c,0,1,-1) +\n       val(c,0,0,-1);\n  m2 = val(c,-1,0,1) + val(c,1,0,1) + val(c,0,-1,1) + val(c,0,1,1) +\n       val(c,0,0,1);\n  m[2][2] = m1 > m2 ? 1. : -1.;\n\n\n  t0 = fabs(m[0][0]) + fabs(m[0][1]) + fabs(m[0][2]);\n  m[0][0] /= t0;\n  m[0][1] /= t0;\n  m[0][2] /= t0;\n\n  t0 = fabs(m[1][0]) + fabs(m[1][1]) + fabs(m[1][2]);\n  m[1][0] /= t0;\n  m[1][1] /= t0;\n  m[1][2] /= t0;\n\n  t0 = fabs(m[2][0]) + fabs(m[2][1]) + fabs(m[2][2]);\n  m[2][0] /= t0;\n  m[2][1] /= t0;\n  m[2][2] /= t0;\n\n\n  t0 = fabs(m[0][0]);\n  t1 = fabs(m[1][1]);\n  t2 = fabs(m[2][2]);\n\n  cn = 0;\n  if (t1 > t0) {\n    t0 = t1;\n    cn = 1;\n  }\n  if (t2 > t0)\n    cn = 2;\n\n\n  m1 = val(c,-1,-1,-1) + val(c,-1,1,-1) + val(c,-1,-1,1) + val(c,-1,1,1) +\n       2.*(val(c,-1,-1,0) + val(c,-1,1,0) + val(c,-1,0,-1) + val(c,-1,0,1)) +\n       4.*val(c,-1,0,0);\n  m2 = val(c,1,-1,-1) + val(c,1,1,-1) + val(c,1,-1,1) + val(c,1,1,1) +\n       2.*(val(c,1,-1,0) + val(c,1,1,0) + val(c,1,0,-1) + val(c,1,0,1)) +\n       4.*val(c,1,0,0);\n  m[3][0] = m1 - m2;\n\n  m1 = val(c,-1,-1,-1) + val(c,-1,-1,1) + val(c,1,-1,-1) + val(c,1,-1,1) +\n       2.*( val(c,-1,-1,0) + val(c,1,-1,0) + val(c,0,-1,-1) + val(c,0,-1,1)) +\n       4.*val(c,0,-1,0);\n  m2 = val(c,-1,1,-1) + val(c,-1,1,1) + val(c,1,1,-1) + val(c,1,1,1) +\n       2.*(val(c,-1,1,0) + val(c,1,1,0) + val(c,0,1,-1) + val(c,0,1,1)) +\n       4.*val(c,0,1,0);\n  m[3][1] = m1 - m2;\n\n  m1 = val(c,-1,-1,-1) + val(c,-1,1,-1) + val(c,1,-1,-1) + val(c,1,1,-1) +\n       2.*(val(c,-1,0,-1) + val(c,1,0,-1) + val(c,0,-1,-1) + val(c,0,1,-1)) +\n       4.*val(c,0,0,-1);\n  m2 = val(c,-1,-1,1) + val(c,-1,1,1) + val(c,1,-1,1) + val(c,1,1,1) +\n       2.*(val(c,-1,0,1) + val(c,1,0,1) + val(c,0,-1,1) + val(c,0,1,1)) +\n       4.*val(c,0,0,1);\n  m[3][2] = m1 - m2;\n\n\n  t0 = fabs(m[3][0]) + fabs(m[3][1]) + fabs(m[3][2]);\n  if (t0 < 1e-30) {\n    coord mxyz = {1., 0., 0.};\n    return mxyz;\n  }\n\n  m[3][0] /= t0;\n  m[3][1] /= t0;\n  m[3][2] /= t0;\n\n\n  t0 = fabs (m[3][0]);\n  t1 = fabs (m[3][1]);\n  t2 = fabs (m[3][2]);\n  if (t1 > t0)\n    t0 = t1;\n  if (t2 > t0)\n    t0 = t2;\n\n  if (fabs(m[cn][cn]) > t0)\n    cn = 3;\n\n\n  coord mxyz = {m[cn][0], m[cn][1], m[cn][2]};\n  return mxyz;\n}","\n    \n// #line 549 \"/home/lisergey/basilisk/src/fractions.h\"\nif (val(c,0,0,0) > 1e-6 && val(c,0,0,0) < 1. - 1e-6) {\n      coord n = mycs (point, c), p;\n      real alpha = plane_alpha (val(c,0,0,0), n);\n      area += pow(Delta, 3 - 1)*plane_area_center (n, alpha, &p);\n    }"})
    {_stencil_val(c,0,0,0); _stencil_val(c,0,0,0); {   
       _stencil_mycs (point, c);     
      _stencil_val(c,0,0,0); 
             
    }        }end_foreach_stencil();
  
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel reduction(+:area)){
#line 548
foreach ()
    if (val(c,0,0,0) > 1e-6 && val(c,0,0,0) < 1. - 1e-6) {
      coord n = mycs (point, c), p;
      double alpha = plane_alpha (val(c,0,0,0), n);
      area += pow(Delta, 3 - 1)*plane_area_center (n, alpha, &p);
    }end_foreach();mpi_all_reduce_array(&area,double,MPI_SUM,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
  
#line 554
{end_tracing("interface_area","/home/lisergey/basilisk/src/fractions.h",554);return area;}
end_tracing("interface_area","/home/lisergey/basilisk/src/fractions.h",555);}
#line 9 "cylinder.c"
#line 1 "navier-stokes/centered.h"
#line 1 "/home/lisergey/basilisk/src/navier-stokes/centered.h"
#line 27 "/home/lisergey/basilisk/src/navier-stokes/centered.h"
#line 1 "./run.h"
#line 1 "/home/lisergey/basilisk/src/run.h"
#line 9 "/home/lisergey/basilisk/src/run.h"
double dt = 1.;

#line 1 "./utils.h"
#line 1 "/home/lisergey/basilisk/src/utils.h"







double DT = 1e30, CFL = 0.5;




struct {

  long nc;

  long tnc;

  double t;

  double speed;

  timer gt;
} perf = {0};





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
#if 1
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
#if 1
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
    "\n# " "Octree"
    ", %d steps, %g CPU, %.4g real, %.3g points.step/s, %d var\n",
    i, s.cpu, s.real, s.speed, (int) (datasize/sizeof(double)));
#if 1
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
  if(!is_constant(cm)){
  
#line 135
foreach_stencil(1,{(NonLocal[]){{"rms","double",(void *)&rms,NULL,0,'+'},{"avg","double",(void *)&avg,NULL,0,'+'},{"volume","double",(void *)&volume,NULL,0,'+'},{"max","double",(void *)&max,NULL,0,'M'},{"cm","scalar",(void *)&cm,NULL,0},{"f","scalar",(void *)&f,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n    \n// #line 137 \"/home/lisergey/basilisk/src/utils.h\"\nif (val(f,0,0,0) != 1e30 && (cube(Delta)*val(cm,0,0,0)) > 0.) {\n      real v = fabs(val(f,0,0,0));\n      if (v > max) max = v;\n      volume += (cube(Delta)*val(cm,0,0,0));\n      avg += (cube(Delta)*val(cm,0,0,0))*v;\n      rms += (cube(Delta)*val(cm,0,0,0))*sq(v);\n    }"}
)
    {_stencil_val(f,0,0,0);_stencil_val(cm,0,0,0); {   
      _stencil_val(f,0,0,0);   
         
_stencil_val(cm,0,0,0); 
       _stencil_val(cm,0,0,0); 
       _stencil_val(cm,0,0,0); 
       
    
#line 143
}       }end_foreach_stencil();
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel  reduction(+:volume)
   reduction(+:rms) reduction(+:avg)reduction(max:max)){
#line 135
foreach(
)
    if (val(f,0,0,0) != 1e30 && (cube(Delta)*val(cm,0,0,0)) > 0.) {
      double v = fabs(val(f,0,0,0));
      if (v > max) max = v;
      volume += (cube(Delta)*val(cm,0,0,0));
      avg += (cube(Delta)*val(cm,0,0,0))*v;
      rms += (cube(Delta)*val(cm,0,0,0))*sq(v);
    }end_foreach();mpi_all_reduce_array(&volume,double,MPI_SUM,1);mpi_all_reduce_array(&rms,double,MPI_SUM,1);mpi_all_reduce_array(&avg,double,MPI_SUM,1);mpi_all_reduce_array(&max,double,MPI_MAX,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 143
}else {double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);
  
#line 135
foreach_stencil(1,{(NonLocal[]){{"rms","double",(void *)&rms,NULL,0,'+'},{"avg","double",(void *)&avg,NULL,0,'+'},{"volume","double",(void *)&volume,NULL,0,'+'},{"max","double",(void *)&max,NULL,0,'M'},{"_const_cm","double",(void *)&_const_cm,NULL,0},{"f","scalar",(void *)&f,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n    \n// #line 137 \"/home/lisergey/basilisk/src/utils.h\"\nif (val(f,0,0,0) != 1e30 && (cube(Delta)*_const_cm) > 0.) {\n      real v = fabs(val(f,0,0,0));\n      if (v > max) max = v;\n      volume += (cube(Delta)*_const_cm);\n      avg += (cube(Delta)*_const_cm)*v;\n      rms += (cube(Delta)*_const_cm)*sq(v);\n    }"}
)
    {_stencil_val(f,0,0,0);; {   
      _stencil_val(f,0,0,0);

;
;
; 
       
    
#line 143
}       }end_foreach_stencil();
  
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel  reduction(+:volume)
   reduction(+:rms) reduction(+:avg)reduction(max:max)){
#line 135
foreach(
)
    if (val(f,0,0,0) != 1e30 && (cube(Delta)*_const_cm) > 0.) {
      double v = fabs(val(f,0,0,0));
      if (v > max) max = v;
      volume += (cube(Delta)*_const_cm);
      avg += (cube(Delta)*_const_cm)*v;
      rms += (cube(Delta)*_const_cm)*sq(v);
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
  if(!is_constant(cm)){
  
#line 163
foreach_stencil(1,{(NonLocal[]){{"min","double",(void *)&min,NULL,0,'m'},{"max","double",(void *)&max,NULL,0,'M'},{"sum2","double",(void *)&sum2,NULL,0,'+'},{"sum","double",(void *)&sum,NULL,0,'+'},{"volume","double",(void *)&volume,NULL,0,'+'},{"f","scalar",(void *)&f,NULL,0},{"cm","scalar",(void *)&cm,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n    \n// #line 165 \"/home/lisergey/basilisk/src/utils.h\"\nif ((cube(Delta)*val(cm,0,0,0)) > 0. && val(f,0,0,0) != 1e30) {\n      volume += (cube(Delta)*val(cm,0,0,0));\n      sum += (cube(Delta)*val(cm,0,0,0))*val(f,0,0,0);\n      sum2 += (cube(Delta)*val(cm,0,0,0))*sq(val(f,0,0,0));\n      if (val(f,0,0,0) > max) max = val(f,0,0,0);\n      if (val(f,0,0,0) < min) min = val(f,0,0,0);\n    }"}
)
    {_stencil_val(cm,0,0,0); _stencil_val(f,0,0,0); {
_stencil_val(cm,0,0,0); 
       _stencil_val(cm,0,0,0);_stencil_val(f,0,0,0); 
       _stencil_val(cm,0,0,0);_stencil_val(f,0,0,0); 
       _stencil_val(f,0,0,0); { _stencil_val(f,0,0,0); }
_stencil_val(f,0,0,0); { _stencil_val(f,0,0,0); }
         
         
    
#line 171
}      }end_foreach_stencil();
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel  reduction(min:min)
   reduction(max:max) reduction(+:volume) reduction(+:sum2)reduction(+:sum)){
#line 163
foreach(
)
    if ((cube(Delta)*val(cm,0,0,0)) > 0. && val(f,0,0,0) != 1e30) {
      volume += (cube(Delta)*val(cm,0,0,0));
      sum += (cube(Delta)*val(cm,0,0,0))*val(f,0,0,0);
      sum2 += (cube(Delta)*val(cm,0,0,0))*sq(val(f,0,0,0));
      if (val(f,0,0,0) > max) max = val(f,0,0,0);
      if (val(f,0,0,0) < min) min = val(f,0,0,0);
    }end_foreach();mpi_all_reduce_array(&min,double,MPI_MIN,1);mpi_all_reduce_array(&max,double,MPI_MAX,1);mpi_all_reduce_array(&volume,double,MPI_SUM,1);mpi_all_reduce_array(&sum2,double,MPI_SUM,1);mpi_all_reduce_array(&sum,double,MPI_SUM,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 171
}else {double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);
  
#line 163
foreach_stencil(1,{(NonLocal[]){{"min","double",(void *)&min,NULL,0,'m'},{"max","double",(void *)&max,NULL,0,'M'},{"sum2","double",(void *)&sum2,NULL,0,'+'},{"sum","double",(void *)&sum,NULL,0,'+'},{"volume","double",(void *)&volume,NULL,0,'+'},{"f","scalar",(void *)&f,NULL,0},{"_const_cm","double",(void *)&_const_cm,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n    \n// #line 165 \"/home/lisergey/basilisk/src/utils.h\"\nif ((cube(Delta)*_const_cm) > 0. && val(f,0,0,0) != 1e30) {\n      volume += (cube(Delta)*_const_cm);\n      sum += (cube(Delta)*_const_cm)*val(f,0,0,0);\n      sum2 += (cube(Delta)*_const_cm)*sq(val(f,0,0,0));\n      if (val(f,0,0,0) > max) max = val(f,0,0,0);\n      if (val(f,0,0,0) < min) min = val(f,0,0,0);\n    }"}
)
    {; _stencil_val(f,0,0,0); {
;
;_stencil_val(f,0,0,0);
;_stencil_val(f,0,0,0); 
       _stencil_val(f,0,0,0); { _stencil_val(f,0,0,0); }
_stencil_val(f,0,0,0); { _stencil_val(f,0,0,0); }
         
         
    
#line 171
}      }end_foreach_stencil();
  
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel  reduction(min:min)
   reduction(max:max) reduction(+:volume) reduction(+:sum2)reduction(+:sum)){
#line 163
foreach(
)
    if ((cube(Delta)*_const_cm) > 0. && val(f,0,0,0) != 1e30) {
      volume += (cube(Delta)*_const_cm);
      sum += (cube(Delta)*_const_cm)*val(f,0,0,0);
      sum2 += (cube(Delta)*_const_cm)*sq(val(f,0,0,0));
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
#line 187 "/home/lisergey/basilisk/src/utils.h"
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
#line 213 "/home/lisergey/basilisk/src/utils.h"
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
#line 237 "/home/lisergey/basilisk/src/utils.h"
void gradients (scalar * f, vector * g)
{
  if (!(list_len(f) == vectors_len(g))) qassert ("/home/lisergey/basilisk/src/utils.h", 239, "list_len(f) == vectors_len(g)");
  foreach_stencil(1,{(NonLocal[]){{"g","vector",(void *)g,NULL,1},{"f","scalar",(void *)f,NULL,1},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},""," \n// #line 240 \"/home/lisergey/basilisk/src/utils.h\"\n{\n    scalar s; vector v;\n    {forin2 (s,v , f,g) {\n      if (s.gradient)\n { {\n\n\n\n\n\n     val_out_(v.x,0,0,0) = s.gradient (val(s,-1,0,0), val(s,0,0,0), val(s,1,0,0))/Delta;\n } \n// #line 244\n{\n\n\n\n\n\n     val_out_(v.y,0,0,0) = s.gradient (val(s,0,-1,0), val(s,0,0,0), val(s,0,1,0))/Delta;\n } \n// #line 244\n{\n\n\n\n\n\n     val_out_(v.z,0,0,0) = s.gradient (val(s,0,0,-1), val(s,0,0,0), val(s,0,0,1))/Delta;\n }}\n      else\n { {\n\n\n\n\n\n     val_out_(v.x,0,0,0) = (val(s,1,0,0) - val(s,-1,0,0))/(2.*Delta);\n } \n// #line 253\n{\n\n\n\n\n\n     val_out_(v.y,0,0,0) = (val(s,0,1,0) - val(s,0,-1,0))/(2.*Delta);\n } \n// #line 253\n{\n\n\n\n\n\n     val_out_(v.z,0,0,0) = (val(s,0,0,1) - val(s,0,0,-1))/(2.*Delta);\n }}\n    } endforin2()}\n  }"}) {
    scalar s; vector v;
    {vector*_i0=g;scalar*_i1= f;if(_i0)for(v=*_i0,s=*_i1;_i0->x.i>= 0;v=*++_i0,s=*++_i1){ {
      if (_attribute[s.i].gradient)
 { {





_stencil_val(s,-1,0,0); _stencil_val(s,0,0,0); _stencil_val(s,1,0,0);





     
#line 250
_stencil_val_a(v.x,0,0,0);   
 } 
#line 244
{





_stencil_val(s,0,-1,0); _stencil_val(s,0,0,0); _stencil_val(s,0,1,0);





     
#line 250
_stencil_val_a(v.y,0,0,0);   
 } 
#line 244
{





_stencil_val(s,0,0,-1); _stencil_val(s,0,0,0); _stencil_val(s,0,0,1);





     
#line 250
_stencil_val_a(v.z,0,0,0);   
 }}
      else
 { {





_stencil_val(s,1,0,0); _stencil_val(s,-1,0,0);





     
#line 259
_stencil_val_a(v.x,0,0,0);   
 } 
#line 253
{





_stencil_val(s,0,1,0); _stencil_val(s,0,-1,0);





     
#line 259
_stencil_val_a(v.y,0,0,0);   
 } 
#line 253
{





_stencil_val(s,0,0,1); _stencil_val(s,0,0,-1);





     
#line 259
_stencil_val_a(v.z,0,0,0);   
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





     val(v.x,0,0,0) = _attribute[s.i].gradient (val(s,-1,0,0), val(s,0,0,0), val(s,1,0,0))/Delta;
 } 
#line 244
{





     val(v.y,0,0,0) = _attribute[s.i].gradient (val(s,0,-1,0), val(s,0,0,0), val(s,0,1,0))/Delta;
 } 
#line 244
{





     val(v.z,0,0,0) = _attribute[s.i].gradient (val(s,0,0,-1), val(s,0,0,0), val(s,0,0,1))/Delta;
 }}
      else
 { {





     val(v.x,0,0,0) = (val(s,1,0,0) - val(s,-1,0,0))/(2.*Delta);
 } 
#line 253
{





     val(v.y,0,0,0) = (val(s,0,1,0) - val(s,0,-1,0))/(2.*Delta);
 } 
#line 253
{





     val(v.z,0,0,0) = (val(s,0,0,1) - val(s,0,0,-1))/(2.*Delta);
 }}
    }}}
  }end_foreach();}
}
#line 280 "/home/lisergey/basilisk/src/utils.h"
void vorticity (const vector u, scalar omega)
{
  if(!is_constant(fm.x) && !is_constant(cm)){
  
#line 282
foreach_stencil(1,{(NonLocal[]){{"cm","scalar",(void *)&cm,NULL,0},{"u","vector",(void *)&u,NULL,0},{"fm","vector",(void *)&fm,NULL,0},{"omega","scalar",(void *)&omega,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n    \n// #line 283 \"/home/lisergey/basilisk/src/utils.h\"\nval_out_(omega,0,0,0) = ((val(fm.x,1,0,0) - val(fm.x,0,0,0))*val(u.y,0,0,0) +\n        val(fm.x,1,0,0)*val(u.y,1,0,0) - val(fm.x,0,0,0)*val(u.y,-1,0,0) -\n        (val(fm.y,0,1,0) - val(fm.y,0,0,0))*val(u.x,0,0,0) +\n        val(fm.y,0,0,0)*val(u.x,0,-1,0) - val(fm.y,0,1,0)*val(u.x,0,1,0))/(2.*(val(cm,0,0,0) + 0.)*Delta);"})
    {_stencil_val(fm.x,1,0,0); _stencil_val(fm.x,0,0,0);_stencil_val(u.y,0,0,0);
        _stencil_val(fm.x,1,0,0);_stencil_val(u.y,1,0,0); _stencil_val(fm.x,0,0,0);_stencil_val(u.y,-1,0,0);
_stencil_val(fm.y,0,1,0); _stencil_val(fm.y,0,0,0);_stencil_val(u.x,0,0,0);
        _stencil_val(fm.y,0,0,0);_stencil_val(u.x,0,-1,0); _stencil_val(fm.y,0,1,0);_stencil_val(u.x,0,1,0);_stencil_val(cm,0,0,0);
#line 283
_stencil_val_a(omega,0,0,0);      
             

}end_foreach_stencil();{
#line 282
foreach()
    val(omega,0,0,0) = ((val(fm.x,1,0,0) - val(fm.x,0,0,0))*val(u.y,0,0,0) +
        val(fm.x,1,0,0)*val(u.y,1,0,0) - val(fm.x,0,0,0)*val(u.y,-1,0,0) -
        (val(fm.y,0,1,0) - val(fm.y,0,0,0))*val(u.x,0,0,0) +
        val(fm.y,0,0,0)*val(u.x,0,-1,0) - val(fm.y,0,1,0)*val(u.x,0,1,0))/(2.*(val(cm,0,0,0) + 0.)*Delta);end_foreach();}}else if(is_constant(fm.x) && !is_constant(cm)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);
  
#line 282
foreach_stencil(1,{(NonLocal[]){{"cm","scalar",(void *)&cm,NULL,0},{"u","vector",(void *)&u,NULL,0},{"_const_fm","_coord",(void *)&_const_fm,NULL,0},{"omega","scalar",(void *)&omega,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n    \n// #line 283 \"/home/lisergey/basilisk/src/utils.h\"\nval_out_(omega,0,0,0) = ((_const_fm.x - _const_fm.x)*val(u.y,0,0,0) +\n        _const_fm.x*val(u.y,1,0,0) - _const_fm.x*val(u.y,-1,0,0) -\n        (_const_fm.y - _const_fm.y)*val(u.x,0,0,0) +\n        _const_fm.y*val(u.x,0,-1,0) - _const_fm.y*val(u.x,0,1,0))/(2.*(val(cm,0,0,0) + 0.)*Delta);"})
    {;;_stencil_val(u.y,0,0,0);
;_stencil_val(u.y,1,0,0);;_stencil_val(u.y,-1,0,0);
;;_stencil_val(u.x,0,0,0);
;_stencil_val(u.x,0,-1,0);;_stencil_val(u.x,0,1,0);_stencil_val(cm,0,0,0);
#line 283
_stencil_val_a(omega,0,0,0);      
             

}end_foreach_stencil();
  {
#line 282
foreach()
    val(omega,0,0,0) = ((_const_fm.x - _const_fm.x)*val(u.y,0,0,0) +
        _const_fm.x*val(u.y,1,0,0) - _const_fm.x*val(u.y,-1,0,0) -
        (_const_fm.y - _const_fm.y)*val(u.x,0,0,0) +
        _const_fm.y*val(u.x,0,-1,0) - _const_fm.y*val(u.x,0,1,0))/(2.*(val(cm,0,0,0) + 0.)*Delta);end_foreach();}}else if(!is_constant(fm.x) && is_constant(cm)){double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);
  
#line 282
foreach_stencil(1,{(NonLocal[]){{"_const_cm","double",(void *)&_const_cm,NULL,0},{"u","vector",(void *)&u,NULL,0},{"fm","vector",(void *)&fm,NULL,0},{"omega","scalar",(void *)&omega,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n    \n// #line 283 \"/home/lisergey/basilisk/src/utils.h\"\nval_out_(omega,0,0,0) = ((val(fm.x,1,0,0) - val(fm.x,0,0,0))*val(u.y,0,0,0) +\n        val(fm.x,1,0,0)*val(u.y,1,0,0) - val(fm.x,0,0,0)*val(u.y,-1,0,0) -\n        (val(fm.y,0,1,0) - val(fm.y,0,0,0))*val(u.x,0,0,0) +\n        val(fm.y,0,0,0)*val(u.x,0,-1,0) - val(fm.y,0,1,0)*val(u.x,0,1,0))/(2.*(_const_cm + 0.)*Delta);"})
    {_stencil_val(fm.x,1,0,0); _stencil_val(fm.x,0,0,0);_stencil_val(u.y,0,0,0);
        _stencil_val(fm.x,1,0,0);_stencil_val(u.y,1,0,0); _stencil_val(fm.x,0,0,0);_stencil_val(u.y,-1,0,0);
_stencil_val(fm.y,0,1,0); _stencil_val(fm.y,0,0,0);_stencil_val(u.x,0,0,0);
        _stencil_val(fm.y,0,0,0);_stencil_val(u.x,0,-1,0); _stencil_val(fm.y,0,1,0);_stencil_val(u.x,0,1,0);;
#line 283
_stencil_val_a(omega,0,0,0);      
             

}end_foreach_stencil();
  {
#line 282
foreach()
    val(omega,0,0,0) = ((val(fm.x,1,0,0) - val(fm.x,0,0,0))*val(u.y,0,0,0) +
        val(fm.x,1,0,0)*val(u.y,1,0,0) - val(fm.x,0,0,0)*val(u.y,-1,0,0) -
        (val(fm.y,0,1,0) - val(fm.y,0,0,0))*val(u.x,0,0,0) +
        val(fm.y,0,0,0)*val(u.x,0,-1,0) - val(fm.y,0,1,0)*val(u.x,0,1,0))/(2.*(_const_cm + 0.)*Delta);end_foreach();}}else {_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);
  
#line 282
foreach_stencil(1,{(NonLocal[]){{"_const_cm","double",(void *)&_const_cm,NULL,0},{"u","vector",(void *)&u,NULL,0},{"_const_fm","_coord",(void *)&_const_fm,NULL,0},{"omega","scalar",(void *)&omega,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n    \n// #line 283 \"/home/lisergey/basilisk/src/utils.h\"\nval_out_(omega,0,0,0) = ((_const_fm.x - _const_fm.x)*val(u.y,0,0,0) +\n        _const_fm.x*val(u.y,1,0,0) - _const_fm.x*val(u.y,-1,0,0) -\n        (_const_fm.y - _const_fm.y)*val(u.x,0,0,0) +\n        _const_fm.y*val(u.x,0,-1,0) - _const_fm.y*val(u.x,0,1,0))/(2.*(_const_cm + 0.)*Delta);"})
    {;;_stencil_val(u.y,0,0,0);
;_stencil_val(u.y,1,0,0);;_stencil_val(u.y,-1,0,0);
;;_stencil_val(u.x,0,0,0);
;_stencil_val(u.x,0,-1,0);;_stencil_val(u.x,0,1,0);;
#line 283
_stencil_val_a(omega,0,0,0);      
             

}end_foreach_stencil();
  {
#line 282
foreach()
    val(omega,0,0,0) = ((_const_fm.x - _const_fm.x)*val(u.y,0,0,0) +
        _const_fm.x*val(u.y,1,0,0) - _const_fm.x*val(u.y,-1,0,0) -
        (_const_fm.y - _const_fm.y)*val(u.x,0,0,0) +
        _const_fm.y*val(u.x,0,-1,0) - _const_fm.y*val(u.x,0,1,0))/(2.*(_const_cm + 0.)*Delta);end_foreach();}}
}





double change (scalar s, scalar sn)
{
  double max = 0.;
  if(!is_constant(cm)){
  
#line 296
foreach_stencil(1,{(NonLocal[]){{"max","double",(void *)&max,NULL,0,'M'},{"sn","scalar",(void *)&sn,NULL,0},{"s","scalar",(void *)&s,NULL,0},{"cm","scalar",(void *)&cm,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},""," \n// #line 296 \"/home/lisergey/basilisk/src/utils.h\"\n{\n    if ((cube(Delta)*val(cm,0,0,0)) > 0.) {\n      real ds = fabs (val(s,0,0,0) - val(sn,0,0,0));\n      if (ds > max)\n max = ds;\n    }\n    val_out_(sn,0,0,0) = val(s,0,0,0);\n  }"}) {
_stencil_val(cm,0,0,0); {     
       _stencil_val(sn,0,0,0);_stencil_val(s,0,0,0);   
       
  
    } 
_stencil_val(s,0,0,0);
       
    
#line 302
_stencil_val_a(sn,0,0,0); 
  }end_foreach_stencil();
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel reduction(max:max)){
#line 296
foreach() {
    if ((cube(Delta)*val(cm,0,0,0)) > 0.) {
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
  
#line 296
foreach_stencil(1,{(NonLocal[]){{"max","double",(void *)&max,NULL,0,'M'},{"sn","scalar",(void *)&sn,NULL,0},{"s","scalar",(void *)&s,NULL,0},{"_const_cm","double",(void *)&_const_cm,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},""," \n// #line 296 \"/home/lisergey/basilisk/src/utils.h\"\n{\n    if ((cube(Delta)*_const_cm) > 0.) {\n      real ds = fabs (val(s,0,0,0) - val(sn,0,0,0));\n      if (ds > max)\n max = ds;\n    }\n    val_out_(sn,0,0,0) = val(s,0,0,0);\n  }"}) {
; {     
       _stencil_val(sn,0,0,0);_stencil_val(s,0,0,0);   
       
  
    } 
_stencil_val(s,0,0,0);
       
    
#line 302
_stencil_val_a(sn,0,0,0); 
  }end_foreach_stencil();
  
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel reduction(max:max)){
#line 296
foreach() {
    if ((cube(Delta)*_const_cm) > 0.) {
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
#line 340 "/home/lisergey/basilisk/src/utils.h"
#define foreach_segment(_S,_p) {\
  coord t = {(_S)[1].x - (_S)[0].x, (_S)[1].y - (_S)[0].y};\
  double norm = sqrt(sq(t.x) + sq(t.y));\
  if (!(norm > 0.)) qassert ("/home/lisergey/basilisk/src/utils.h", 343, "norm > 0.");\
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

#line 410 "/home/lisergey/basilisk/src/utils.h"
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
#line 1 "/home/lisergey/basilisk/src/output.h"
#line 37 "/home/lisergey/basilisk/src/output.h"
     
void output_field (scalar * list,
     FILE * fp,
     int n,
     bool linear,
     coord box[2])
{tracing("output_field","/home/lisergey/basilisk/src/output.h",38);
  n++;
  int len = list_len (list);
  double Delta = 0.999999*(box[1].x - box[0].x)/(n - 1);
  int ny = (box[1].y - box[0].y)/Delta + 1;
  double ** field = (double **) matrix_new (n, ny, len*sizeof(double)), * v = field[0];
  for (int i = 0; i < n*ny*len; i++, v++)
    *v = 1e30;
  coord box1[2] = {{box[0].x - Delta/2., box[0].y - Delta/2.},
     {box[0].x + (n - 0.5)*Delta, box[0].y + (ny - 0.5)*Delta}};
  coord cn = {n, ny}, p;

  v = field[0];
  foreach_region_stencil (1,{(NonLocal[]){{"linear","bool",(void *)&linear,NULL,0},{"list","scalar",(void *)list,NULL,1},{"field","double",(void *)field,NULL,2},{"len","int",(void *)&len,NULL,0},{"ny","int",(void *)&ny,NULL,0},{"n","int",(void *)&n,NULL,0},{"cn","coord",(void *)&cn,NULL,0},{"box1","coord",(void *)box1,(int[]){2,0},0},{"p","coord",(void *)&p,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"\n\n\n// #line 805 \"/home/lisergey/basilisk/src/grid/cartesian-common.h\"\nstatic real interpolate_linear (Point point, scalar v,\n      real xp, real yp, real zp)\n{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;\n// # 823 \"/home/lisergey/basilisk/src/grid/cartesian-common.h\"\n  x = (xp - x)/Delta - v.d.x/2.;\n  y = (yp - y)/Delta - v.d.y/2.;\n  z = (zp - z)/Delta - v.d.z/2.;\n  int i = sign(x), j = sign(y), k = sign(z);\n  x = fabs(x); y = fabs(y); z = fabs(z);\n\n  return (((val(v,0,0,0)*(1. - x) + val(v,i,0,0)*x)*(1. - y) +\n    (val(v,0,j,0)*(1. - x) + val(v,i,j,0)*x)*y)*(1. - z) +\n   ((val(v,0,0,k)*(1. - x) + val(v,i,0,k)*x)*(1. - y) +\n    (val(v,0,j,k)*(1. - x) + val(v,i,j,k)*x)*y)*z);\n\n}","\n\n\n\n  \n// #line 60 \"/home/lisergey/basilisk/src/output.h\"\n{\n    real ** alias = field;\n    int i = (p.x - box1[0].x)/(box1[1].x - box1[0].x)*cn.x;\n    int j = (p.y - box1[0].y)/(box1[1].y - box1[0].y)*cn.y;\n    int k = 0;\n    {forin (scalar, s , list)\n      alias[i][len*j + k++] = linear ? interpolate_linear (point, s, p.x, p.y, p.z) : val(s,0,0,0); endforin()}\n  }"})



  {                     
    
    
    
    
    {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      { _stencil_interpolate_linear (point, s, NULL, NULL, NULL); _stencil_val(s,0,0,0);      }}}
  }end_foreach_region_stencil();
  
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel  reduction(min:v[:n*ny*len])){
#line 56
foreach_region (p, box1, cn)



  {
    double ** alias = field;
    int i = (p.x - box1[0].x)/(box1[1].x - box1[0].x)*cn.x;
    int j = (p.y - box1[0].y)/(box1[1].y - box1[0].y)*cn.y;
    int k = 0;
    {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      alias[i][len*j + k++] = linear ? interpolate_linear (point, s, p.x, p.y, p.z) : val(s,0,0,0);}}
  }end_foreach_region();mpi_all_reduce_array(v,double,MPI_MIN,n*ny*len);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}

  
#line 69
if (pid() == 0) {
    fprintf (fp, "# 1:x 2:y");
    int i = 3;
    {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      fprintf (fp, " %d:%s", i++, _attribute[s.i].name);}}
    fputc('\n', fp);
    for (int i = 0; i < n; i++) {
      double x = Delta*i + box[0].x;
      for (int j = 0; j < ny; j++) {
 double y = Delta*j + box[0].y;

 fprintf (fp, "%g %g", x, y);
 int k = 0;
 {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
   fprintf (fp, " %g", field[i][len*j + k++]);}}
 fputc ('\n', fp);
      }
      fputc ('\n', fp);
    }
    fflush (fp);
  }

  matrix_free (field);
end_tracing("output_field","/home/lisergey/basilisk/src/output.h",92);}
#line 120 "/home/lisergey/basilisk/src/output.h"
     
void output_matrix (scalar f, FILE * fp, int n, bool linear)
{tracing("output_matrix","/home/lisergey/basilisk/src/output.h",121);
  float fn = n;
  float Delta = (float) L0/fn;
  fwrite (&fn, sizeof(float), 1, fp);
  for (int j = 0; j < n; j++) {
    float yp = (float) (Delta*j + X0 + Delta/2.);
    fwrite (&yp, sizeof(float), 1, fp);
  }
  for (int i = 0; i < n; i++) {
    float xp = (float) (Delta*i + X0 + Delta/2.);
    fwrite (&xp, sizeof(float), 1, fp);
    for (int j = 0; j < n; j++) {
      float yp = (float)(Delta*j + Y0 + Delta/2.), v;
      v = interpolate (f, xp, yp
#line 837 "/home/lisergey/basilisk/src/grid/cartesian-common.h"
, 0.
#line 135 "/home/lisergey/basilisk/src/output.h"
, linear);
      fwrite (&v, sizeof(float), 1, fp);
    }
  }
  fflush (fp);
end_tracing("output_matrix","/home/lisergey/basilisk/src/output.h",140);}
#line 149 "/home/lisergey/basilisk/src/output.h"
typedef void (* Colormap) (double cmap[127][3]);

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
} Color;

Color colormap_color (double cmap[127][3],
        double val, double min, double max)
{
  Color c;
  if (val == 1e30) {
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
  if (!(i >= 0 && i < 127 - 1)) qassert ("/home/lisergey/basilisk/src/output.h", 281, "i >= 0 && i < NCMAP - 1");
  unsigned char * c1 = (unsigned char *) &c;
  for (int j = 0; j < 3; j++)
    c1[j] = 255*(cmap[i][j]*(1. - coef) + cmap[i + 1][j]*coef);
  return c;
}
#line 300 "/home/lisergey/basilisk/src/output.h"
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
    pfree (open_image_data.names[i],__func__,__FILE__,__LINE__);
  }
  pfree (open_image_data.fp,__func__,__FILE__,__LINE__);
  pfree (open_image_data.names,__func__,__FILE__,__LINE__);
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

    MPI_Abort (MPI_COMM_WORLD, 1);

    exit (1);
  }
  return fp;
}

FILE * open_image (const char * file, const char * options)
{
  if (!(pid() == 0)) qassert ("/home/lisergey/basilisk/src/output.h", 382, "pid() == 0");
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
    open_image_data.names = (char * *) prealloc (open_image_data.names, (open_image_data.n)*sizeof(char *),__func__,__FILE__,__LINE__);
    open_image_data.names[open_image_data.n - 1] = pstrdup (file,__func__,__FILE__,__LINE__);

    if (options) {
      strcat (command, " ");
      strcat (command, options);
    }
    strcat (command, !strcmp (ext, ".mp4") ? " " : " > ");
    strcat (command, file);
    open_image_data.fp = (FILE * *) prealloc (open_image_data.fp, (open_image_data.n)*sizeof(FILE *),__func__,__FILE__,__LINE__);
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
  if (!(pid() == 0)) qassert ("/home/lisergey/basilisk/src/output.h", 456, "pid() == 0");
  if (is_animation (file)) {
    if (!open_image_lookup (file))
      fclose (fp);
  }
  else if (which ("convert"))
    qpclose (fp);
  else
    fclose (fp);
}
#line 534 "/home/lisergey/basilisk/src/output.h"
     
void output_ppm (scalar f,
   FILE * fp,
   int n,
   char * file,
   double min, double max, double spread,
   double z,
   bool linear,
   coord box[2],
   scalar mask,
   Colormap map,
   char * opt)
{tracing("output_ppm","/home/lisergey/basilisk/src/output.h",535);

  if (!min && !max) {
    stats s = statsf (f);
    if (spread < 0.)
      min = s.min, max = s.max;
    else {
      double avg = s.sum/s.volume;
      min = avg - spread*s.stddev; max = avg + spread*s.stddev;
    }
  }
  box[0].z = z, box[1].z = z;

  coord cn = {n}, p;
  double delta = (box[1].x - box[0].x)/n;
  cn.y = (int)((box[1].y - box[0].y)/delta);
  if (((int)cn.y) % 2) cn.y++;

  Color ** ppm = (Color **) matrix_new (cn.y, cn.x, sizeof(Color));
  unsigned char * ppm0 = &ppm[0][0].r;
  int len = 3*cn.x*cn.y;
  memset (ppm0, 0, len*sizeof (unsigned char));
  double cmap[127][3];
  (* map) (cmap);


  foreach_region_stencil (1,{(NonLocal[]){{"max","double",(void *)&max,NULL,0},{"min","double",(void *)&min,NULL,0},{"cmap","double",(void *)cmap,(int[]){1273,0},0},{"ferr","not implemented yet",(void *)ferr,NULL,1},{"ppm","Color",(void *)ppm,NULL,2},{"f","scalar",(void *)&f,NULL,0},{"linear","bool",(void *)&linear,NULL,0},{"mask","scalar",(void *)&mask,NULL,0},{"len","int",(void *)&len,NULL,0},{"cn","coord",(void *)&cn,NULL,0},{"box","coord",(void *)box,(int[]){2,0},0},{"p","coord",(void *)&p,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"\n@ def not_mpi_compatible()\ndo {\n  if (npe() > 1) {\n    fprintf (ferr, \"%s() is not compatible with MPI (yet)\\n\", __func__);\n    exit (1);\n  }\n} while(0)\n@\n@ define system(command) (pid() == 0 ? system(command) : 0)\n@else\n@ define qstderr() stderr\n@ define qstdout() stdout\n@ define ferr stderr\n@ define fout stdout\n@ define not_mpi_compatible()\n@endif\n\n\n\n\n// #line 93 \"/home/lisergey/basilisk/src/common.h\"\nstatic inline void qassert (const char * file, int line, const char * cond) {\n  fprintf (ferr, \"%s:%d: Assertion `%s' failed.\\n\", file, line, cond);\n  abort();\n}\n\n\n// #line 261 \"/home/lisergey/basilisk/src/output.h\"\nColor colormap_color (real cmap[127][3],\n        real val, real min, real max)\n{\n  Color c;\n  if (val == 1e30) {\n    c.r = c.g = c.b = 0;\n    return c;\n  }\n  int i;\n  real coef;\n  if (max != min)\n    val = (val - min)/(max - min);\n  else\n    val = 0.;\n  if (val <= 0.) i = 0, coef = 0.;\n  else if (val >= 1.) i = 127 - 2, coef = 1.;\n  else {\n    i = val*(127 - 1);\n    coef = val*(127 - 1) - i;\n  }\n  if (!(i >= 0 && i < 127 - 1)) qassert (\"/home/lisergey/basilisk/src/output.h\", 281, \"i >= 0 && i < NCMAP - 1\");\n  unsigned char * c1 = (unsigned char *) &c;\n  for (int j = 0; j < 3; j++)\n    c1[j] = 255*(cmap[i][j]*(1. - coef) + cmap[i + 1][j]*coef);\n  return c;\n}\n\n\n// #line 805 \"/home/lisergey/basilisk/src/grid/cartesian-common.h\"\nstatic real interpolate_linear (Point point, scalar v,\n      real xp, real yp, real zp)\n{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;\n// # 823 \"/home/lisergey/basilisk/src/grid/cartesian-common.h\"\n  x = (xp - x)/Delta - v.d.x/2.;\n  y = (yp - y)/Delta - v.d.y/2.;\n  z = (zp - z)/Delta - v.d.z/2.;\n  int i = sign(x), j = sign(y), k = sign(z);\n  x = fabs(x); y = fabs(y); z = fabs(z);\n\n  return (((val(v,0,0,0)*(1. - x) + val(v,i,0,0)*x)*(1. - y) +\n    (val(v,0,j,0)*(1. - x) + val(v,i,j,0)*x)*y)*(1. - z) +\n   ((val(v,0,0,k)*(1. - x) + val(v,i,0,k)*x)*(1. - y) +\n    (val(v,0,j,k)*(1. - x) + val(v,i,j,k)*x)*y)*z);\n\n}","\n\n\n\n  \n// #line 576 \"/home/lisergey/basilisk/src/output.h\"\n{\n    real v;\n    if (mask.i >= 0) {\n      if (linear) {\n real m = interpolate_linear (point, mask, p.x, p.y, p.z);\n if (m < 0.)\n   v = 1e30;\n else\n   v = interpolate_linear (point, f, p.x, p.y, p.z);\n      }\n      else {\n if (val(mask,0,0,0) < 0.)\n   v = 1e30;\n else\n   v = val(f,0,0,0);\n      }\n    }\n    else if (linear)\n      v = interpolate_linear (point, f, p.x, p.y, p.z);\n    else\n      v = val(f,0,0,0);\n    int i = (p.x - box[0].x)/(box[1].x - box[0].x)*cn.x;\n    int j = (p.y - box[0].y)/(box[1].y - box[0].y)*cn.y;\n    Color ** alias = ppm;\n    alias[(int)cn.y - 1 - j][i] = colormap_color (cmap, v, min, max);\n  }"})



  { 
    
    if (mask.i >= 0) {
      if (linear) {  
  _stencil_interpolate_linear (point, mask, NULL, NULL, NULL);
{ 
    
   
{ _stencil_interpolate_linear (point, f, NULL, NULL, NULL); }}
    
 
      
#line 585
}
      else {
_stencil_val(mask,0,0,0);{
     
   
{ _stencil_val(f,0,0,0); }}
    
 
      
#line 591
}
    }
    else if (linear)
      { _stencil_interpolate_linear (point, f, NULL, NULL, NULL); }
    else
      { _stencil_val(f,0,0,0); }                  
    
    
         
         
  }end_foreach_region_stencil();


  
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel  reduction(max:ppm0[:len])){
#line 572
foreach_region (p, box, cn)



  {
    double v;
    if (mask.i >= 0) {
      if (linear) {
 double m = interpolate_linear (point, mask, p.x, p.y, p.z);
 if (m < 0.)
   v = 1e30;
 else
   v = interpolate_linear (point, f, p.x, p.y, p.z);
      }
      else {
 if (val(mask,0,0,0) < 0.)
   v = 1e30;
 else
   v = val(f,0,0,0);
      }
    }
    else if (linear)
      v = interpolate_linear (point, f, p.x, p.y, p.z);
    else
      v = val(f,0,0,0);
    int i = (p.x - box[0].x)/(box[1].x - box[0].x)*cn.x;
    int j = (p.y - box[0].y)/(box[1].y - box[0].y)*cn.y;
    Color ** alias = ppm;
    alias[(int)cn.y - 1 - j][i] = colormap_color (cmap, v, min, max);
  }end_foreach_region();mpi_all_reduce_array(ppm0,unsigned char,MPI_MAX,len);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}

  
#line 603
if (pid() == 0) {
    if (file)
      fp = open_image (file, opt);

    fprintf (fp, "P6\n%g %g 255\n", cn.x, cn.y);
    fwrite (ppm0, sizeof(unsigned char), 3*cn.x*cn.y, fp);

    if (file)
      close_image (file, fp);
    else
      fflush (fp);
  }

  matrix_free (ppm);
end_tracing("output_ppm","/home/lisergey/basilisk/src/output.h",617);}
#line 649 "/home/lisergey/basilisk/src/output.h"
     
void output_grd (scalar f,
   FILE * fp,
   double Delta,
   bool linear,
   double box[2][2],
   scalar mask)
{tracing("output_grd","/home/lisergey/basilisk/src/output.h",650);
  int nx = (box[1][0] - box[0][0])/Delta;
  int ny = (box[1][1] - box[0][1])/Delta;


  fprintf (fp, "ncols          %d\n", nx);
  fprintf (fp, "nrows          %d\n", ny);
  fprintf (fp, "xllcorner      %g\n", box[0][0]);
  fprintf (fp, "yllcorner      %g\n", box[0][1]);
  fprintf (fp, "cellsize       %g\n", Delta);
  fprintf (fp, "nodata_value   -9999\n");


  for (int j = ny-1; j >= 0; j--) {
    double yp = Delta*j + box[0][1] + Delta/2.;
    for (int i = 0; i < nx; i++) {
      double xp = Delta*i + box[0][0] + Delta/2., v;
      if (mask.i >= 0) {
 double m = interpolate (mask, xp, yp
#line 837 "/home/lisergey/basilisk/src/grid/cartesian-common.h"
, 0.
#line 674 "/home/lisergey/basilisk/src/output.h"
, linear);
 if (m < 0.)
   v = 1e30;
 else
   v = interpolate (f, xp, yp
#line 837 "/home/lisergey/basilisk/src/grid/cartesian-common.h"
, 0.
#line 678 "/home/lisergey/basilisk/src/output.h"
, linear);
      }
      else
 v = interpolate (f, xp, yp
#line 837 "/home/lisergey/basilisk/src/grid/cartesian-common.h"
, 0.
#line 681 "/home/lisergey/basilisk/src/output.h"
, linear);
      if (v == 1e30)
 fprintf (fp, "-9999 ");
      else
 fprintf (fp, "%f ", v);
    }
    fprintf (fp, "\n");
  }

  fflush (fp);
end_tracing("output_grd","/home/lisergey/basilisk/src/output.h",691);}
#line 718 "/home/lisergey/basilisk/src/output.h"
static char * replace (const char * input, int target, int with,
         bool translate)
{
  if (translate) {
    if (!strcmp (input, "u.x"))
      return pstrdup ("U",__func__,__FILE__,__LINE__);
    if (!strcmp (input, "u.y"))
      return pstrdup ("V",__func__,__FILE__,__LINE__);
    if (!strcmp (input, "u.z"))
      return pstrdup ("W",__func__,__FILE__,__LINE__);
  }
  char * name = pstrdup (input,__func__,__FILE__,__LINE__), * i = name;
  while (*i != '\0') {
    if (*i == target)
      *i = with;
    i++;
  }
  return name;
}

     
void output_gfs (FILE * fp,
   scalar * list,
   char * file,
   bool translate)
{tracing("output_gfs","/home/lisergey/basilisk/src/output.h",739);
  char * fname = file;

#if 1



  FILE * sfp = fp;
  if (file == NULL) {
    long pid = getpid();
    MPI_Bcast (&pid, 1, MPI_LONG, 0, MPI_COMM_WORLD);
    fname = ((char *) pmalloc ((80)*sizeof(char),__func__,__FILE__,__LINE__));
    snprintf (fname, 80, ".output-%ld", pid);
    fp = NULL;
  }
#endif

  bool opened = false;
  if (fp == NULL) {
    if (fname == NULL)
      fp = fout;
    else if (!(fp = fopen (fname, "w"))) {
      perror (fname);
      exit (1);
    }
    else
      opened = true;
  }

  scalar * slist = list ? list : list_copy (all);

  restriction (slist);
  fprintf (fp,
    "1 0 GfsSimulation GfsBox GfsGEdge { binary = 1"
    " x = %g y = %g ",
    0.5 + X0/L0, 0.5 + Y0/L0);

  fprintf (fp, "z = %g ", 0.5 + Z0/L0);


  if (slist != NULL && slist[0].i != -1) {
    scalar s = slist[0];
    char * name = replace (_attribute[s.i].name, '.', '_', translate);
    fprintf (fp, "variables = %s", name);
    pfree (name,__func__,__FILE__,__LINE__);
    for (int i = 1; i < list_len(slist); i++) {
      scalar s = slist[i];
      if (_attribute[s.i].name) {
 char * name = replace (_attribute[s.i].name, '.', '_', translate);
 fprintf (fp, ",%s", name);
 pfree (name,__func__,__FILE__,__LINE__);
      }
    }
    fprintf (fp, " ");
  }
  fprintf (fp, "} {\n");
  fprintf (fp, "  Time { t = %g }\n", t);
  if (L0 != 1.)
    fprintf (fp, "  PhysicalParams { L = %g }\n", L0);
  fprintf (fp, "  VariableTracerVOF f\n");
  fprintf (fp, "}\nGfsBox { x = 0 y = 0 z = 0 } {\n");

#if 1
  long header;
  if ((header = ftell (fp)) < 0) {
    perror ("output_gfs(): error in header");
    exit (1);
  }
  int cell_size = sizeof(unsigned) + sizeof(double);
  {scalar*_i=(scalar*)( slist);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    if (_attribute[s.i].name)
      cell_size += sizeof(double);}}
  scalar index = new_scalar("index");
  size_t total_size = header + (z_indexing (index, false) + 1)*cell_size;
#endif



  {foreach_cell() {
#if 1
    if (is_local(cell))
#endif
    {
#if 1
      if (fseek (fp, header + val(index,0,0,0)*cell_size, SEEK_SET) < 0) {
 perror ("output_gfs(): error while seeking");
 exit (1);
      }
#endif
      unsigned flags =
 level == 0 ? 0 :
#line 842 "/home/lisergey/basilisk/src/output.h"
      child.x == -1 && child.y == -1 && child.z == -1 ? 0 :
 child.x == -1 && child.y == -1 && child.z == 1 ? 1 :
 child.x == -1 && child.y == 1 && child.z == -1 ? 2 :
 child.x == -1 && child.y == 1 && child.z == 1 ? 3 :
 child.x == 1 && child.y == -1 && child.z == -1 ? 4 :
 child.x == 1 && child.y == -1 && child.z == 1 ? 5 :
 child.x == 1 && child.y == 1 && child.z == -1 ? 6 :
 7;

      if (is_leaf(cell))
 flags |= (1 << 4);
      fwrite (&flags, sizeof (unsigned), 1, fp);
      double a = -1;
      fwrite (&a, sizeof (double), 1, fp);
      {scalar*_i=(scalar*)( slist);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
 if (_attribute[s.i].name) {
   if (_attribute[s.i].v.x.i >= 0) {




     if (_attribute[s.i].v.x.i == s.i) {
       s = _attribute[s.i].v.y;
       a = is_local(cell) && val(s,0,0,0) != 1e30 ? val(s,0,0,0) : (double) DBL_MAX;
     }
     else if (_attribute[s.i].v.y.i == s.i) {
       s = _attribute[s.i].v.x;
       a = is_local(cell) && val(s,0,0,0) != 1e30 ? - val(s,0,0,0) : (double) DBL_MAX;
     }


     else
       a = is_local(cell) && val(s,0,0,0) != 1e30 ? val(s,0,0,0) : (double) DBL_MAX;

   }
   else
     a = is_local(cell) && val(s,0,0,0) != 1e30 ? val(s,0,0,0) : (double) DBL_MAX;
   fwrite (&a, sizeof (double), 1, fp);
 }}}
    }
    if (is_leaf(cell))
      continue;
  }end_foreach_cell();}

#if 1
  delete (((scalar[]){index,{-1}}));
  if (!pid() && fseek (fp, total_size, SEEK_SET) < 0) {
    perror ("output_gfs(): error while finishing");
    exit (1);
  }
  if (!pid())
#endif
    fputs ("}\n", fp);
  fflush (fp);

  if (!list)
    pfree (slist,__func__,__FILE__,__LINE__);
  if (opened)
    fclose (fp);

#if 1
  if (file == NULL) {
    MPI_Barrier (MPI_COMM_WORLD);
    if (pid() == 0) {
      if (sfp == NULL)
 sfp = fout;
      fp = fopen (fname, "r");
      size_t l;
      unsigned char buffer[8192];
      while ((l = fread (buffer, 1, 8192, fp)) > 0)
 fwrite (buffer, 1, l, sfp);
      fflush (sfp);
      remove (fname);
    }
    pfree (fname,__func__,__FILE__,__LINE__);
  }
#endif
end_tracing("output_gfs","/home/lisergey/basilisk/src/output.h",919);}
#line 943 "/home/lisergey/basilisk/src/output.h"
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

#if !1
     
void dump (const char * file,
    scalar * list,
    FILE * fp,
    bool unbuffered)
{tracing("dump","/home/lisergey/basilisk/src/output.h",989);
  char * name = NULL;
  if (!fp) {
    name = (char *) pmalloc (strlen(file) + 2,__func__,__FILE__,__LINE__);
    strcpy (name, file);
    if (!unbuffered)
      strcat (name, "~");
    if ((fp = fopen (name, "w")) == NULL) {
      perror (name);
      exit (1);
    }
  }
  if (!(fp)) qassert ("/home/lisergey/basilisk/src/output.h", 1005, "fp");

  scalar * dlist = dump_list (list);
  scalar  size=new_scalar("size");
  scalar * slist = list_concat (((scalar[]){size,{-1}}), dlist); pfree (dlist,__func__,__FILE__,__LINE__);
  struct DumpHeader header = { t, list_len(slist), iter, depth(), npe(),
          dump_version };
  dump_header (fp, &header, slist);

  subtree_size (size, false);

  {foreach_cell() {
    unsigned flags = is_leaf(cell) ? leaf : 0;
    if (fwrite (&flags, sizeof(unsigned), 1, fp) < 1) {
      perror ("dump(): error while writing flags");
      exit (1);
    }
    {scalar*_i=(scalar*)( slist);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      if (fwrite (&val(s,0,0,0), sizeof(double), 1, fp) < 1) {
 perror ("dump(): error while writing scalars");
 exit (1);
      }}}
    if (is_leaf(cell))
      continue;
  }end_foreach_cell();}

  pfree (slist,__func__,__FILE__,__LINE__);
  if (file) {
    fclose (fp);
    if (!unbuffered)
      rename (name, file);
    pfree (name,__func__,__FILE__,__LINE__);
  }delete((scalar*)((scalar[]){size,{-1}}));
end_tracing("dump","/home/lisergey/basilisk/src/output.h",1038);}
#else
     
void dump (const char * file,
    scalar * list,
    FILE * fp,
    bool unbuffered)
{tracing("dump","/home/lisergey/basilisk/src/output.h",1041);
  if (fp != NULL || file == NULL) {
    fprintf (ferr, "dump(): must specify a file name when using MPI\n");
    exit(1);
  }

  char name[strlen(file) + 2];
  strcpy (name, file);
  if (!unbuffered)
    strcat (name, "~");
  FILE * fh = fopen (name, "w");
  if (fh == NULL) {
    perror (name);
    exit (1);
  }

  scalar * dlist = dump_list (list);
  scalar  size=new_scalar("size");
  scalar * slist = list_concat (((scalar[]){size,{-1}}), dlist); pfree (dlist,__func__,__FILE__,__LINE__);
  struct DumpHeader header = { t, list_len(slist), iter, depth(), npe(),
          dump_version };







  if (pid() == 0)
    dump_header (fh, &header, slist);

  scalar index = {-1};

  index = new_scalar("index");
  z_indexing (index, false);
  int cell_size = sizeof(unsigned) + header.len*sizeof(double);
  int sizeofheader = sizeof(header) + 4*sizeof(double);
  {scalar*_i=(scalar*)( slist);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
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
      {scalar*_i=(scalar*)( slist);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
 fwrite (&val(s,0,0,0), 1, sizeof(double), fh);}}
      pos += cell_size;
    }
    if (is_leaf(cell))
      continue;
  }end_foreach_cell();}

  delete (((scalar[]){index,{-1}}));

  pfree (slist,__func__,__FILE__,__LINE__);
  fclose (fh);
  if (!unbuffered && pid() == 0)
    rename (name, file);delete((scalar*)((scalar[]){size,{-1}}));
end_tracing("dump","/home/lisergey/basilisk/src/output.h",1112);}
#endif

     
bool restore (const char * file,
       scalar * list,
       FILE * fp)
{tracing("restore","/home/lisergey/basilisk/src/output.h",1116);
  if (!fp && (fp = fopen (file, "r")) == NULL)
    {end_tracing("restore","/home/lisergey/basilisk/src/output.h",1121);return false;}
  if (!(fp)) qassert ("/home/lisergey/basilisk/src/output.h", 1122, "fp");

  struct DumpHeader header = {0};
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
#line 1157 "/home/lisergey/basilisk/src/output.h"
  bool restore_all = (list == all);
  scalar * slist = dump_list (list ? list : all);
  if (header.version == 161020) {
    if (header.len - 1 != list_len (slist)) {
      fprintf (ferr,
        "restore(): error: the list lengths don't match: "
        "%ld (file) != %d (code)\n",
        header.len - 1, list_len (slist));
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
 {scalar*_i=(scalar*)( slist);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
   if (!strcmp (_attribute[s.i].name, name)) {
     input = list_append (input, s);
     found = true; break;
   }}}
 if (!found) {
   if (restore_all) {
     scalar s = new_scalar("s");
     pfree (_attribute[s.i].name,__func__,__FILE__,__LINE__);
     _attribute[s.i].name = pstrdup (name,__func__,__FILE__,__LINE__);
     input = list_append (input, s);
   }
   else
     input = list_append (input, (scalar){INT_MAX});
 }
      }
    }
    pfree (slist,__func__,__FILE__,__LINE__);
    slist = input;

    double o[4];
    if (fread (o, sizeof(double), 4, fp) < 4) {
      fprintf (ferr, "restore(): error: expecting coordinates\n");
      exit (1);
    }
    origin (o[0], o[1], o[2]);
    size (o[3]);
  }
#line 1232 "/home/lisergey/basilisk/src/output.h"
  scalar * listm = is_constant(cm) ? NULL : (scalar *)((vector[]){fm,{{-1},{-1},{-1}}});

  restore_mpi (fp, slist);
#line 1262 "/home/lisergey/basilisk/src/output.h"
  scalar * other = NULL;
  {scalar*_i=(scalar*)( all);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    if (!list_lookup (slist, s) && !list_lookup (listm, s))
      other = list_append (other, s);}}
  reset (other, 0.);
  pfree (other,__func__,__FILE__,__LINE__);

  pfree (slist,__func__,__FILE__,__LINE__);
  if (file)
    fclose (fp);


  while (iter < header.i && events (false))
    iter = inext;
  events (false);
  while (t < header.t && events (false))
    t = tnext;
  t = header.t;
  events (false);

  {end_tracing("restore","/home/lisergey/basilisk/src/output.h",1282);return true;}
end_tracing("restore","/home/lisergey/basilisk/src/output.h",1283);}
#line 431 "/home/lisergey/basilisk/src/utils.h"
#line 12 "/home/lisergey/basilisk/src/run.h"

     
void run (void)
{tracing("run","/home/lisergey/basilisk/src/run.h",14);
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
end_tracing("run","/home/lisergey/basilisk/src/run.h",37);}




static int defaults_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i = 0)!=0;*ip=i;*tp=t;return ret;}





#line 42
      static int defaults(const int i,const double t,Event *_ev){tracing("defaults","/home/lisergey/basilisk/src/run.h",42); {
  display ("box();"
#line 1425 "/home/lisergey/basilisk/src/common.h"
, false
#line 43 "/home/lisergey/basilisk/src/run.h"
);
}{end_tracing("defaults","/home/lisergey/basilisk/src/run.h",44);return 0;}end_tracing("defaults","/home/lisergey/basilisk/src/run.h",44);}





static int cleanup_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(t = TEND_EVENT)!=0;*ip=i;*tp=t;return ret;}






#line 50
      static int cleanup(const int i,const double t,Event *_ev){tracing("cleanup","/home/lisergey/basilisk/src/run.h",50); {
  display ("", true);
}{end_tracing("cleanup","/home/lisergey/basilisk/src/run.h",52);return 0;}end_tracing("cleanup","/home/lisergey/basilisk/src/run.h",52);}
#line 28 "/home/lisergey/basilisk/src/navier-stokes/centered.h"
#line 1 "./timestep.h"
#line 1 "/home/lisergey/basilisk/src/timestep.h"

double timestep (const vector u, double dtmax)
{
  static double previous = 0.;
  dtmax /= CFL;
  if(!is_constant(cm)){
  
#line 6
foreach_face_stencil(1,{(NonLocal[]){{"dtmax","double",(void *)&dtmax,NULL,0,'m'},{"cm","scalar",(void *)&cm,NULL,0},{"u","vector",(void *)&u,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n// #line 6 \"/home/lisergey/basilisk/src/timestep.h\"\n{is_face_x(){\n    if (val(u.x,0,0,0) != 0.) {\n      real dt = Delta/fabs(val(u.x,0,0,0));\n\n\n\n\n      dt *= val(cm,0,0,0);\n\n      if (dt < dtmax) dtmax = dt;\n    }}end_is_face_x()\n// #line 6\nis_face_y(){\n    if (val(u.y,0,0,0) != 0.) {\n      real dt = Delta/fabs(val(u.y,0,0,0));\n\n\n\n\n      dt *= val(cm,0,0,0);\n\n      if (dt < dtmax) dtmax = dt;\n    }}end_is_face_y()\n// #line 6\nis_face_z(){\n    if (val(u.z,0,0,0) != 0.) {\n      real dt = Delta/fabs(val(u.z,0,0,0));\n\n\n\n\n      dt *= val(cm,0,0,0);\n\n      if (dt < dtmax) dtmax = dt;\n    }}end_is_face_z()}"}){_stencil_is_face_x(){
    {_stencil_val(u.x,0,0,0); {   
      _stencil_val(u.x,0,0,0); 




_stencil_val(cm,0,0,0);   




       

         
    
#line 16
}   }}end__stencil_is_face_x()
#line 6
_stencil_is_face_y(){
    {_stencil_val(u.y,0,0,0); {   
      _stencil_val(u.y,0,0,0); 




_stencil_val(cm,0,0,0);   




       

         
    
#line 16
}   }}end__stencil_is_face_y()
#line 6
_stencil_is_face_z(){
    {_stencil_val(u.z,0,0,0); {   
      _stencil_val(u.z,0,0,0); 




_stencil_val(cm,0,0,0);   




       

         
    
#line 16
}   }}end__stencil_is_face_z()}end_foreach_face_stencil();
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel reduction(min:dtmax)){
#line 6
foreach_face_generic(){is_face_x(){
    if (val(u.x,0,0,0) != 0.) {
      double dt = Delta/fabs(val(u.x,0,0,0));




      dt *= val(cm,0,0,0);

      if (dt < dtmax) dtmax = dt;
    }}end_is_face_x()
#line 6
is_face_y(){
    if (val(u.y,0,0,0) != 0.) {
      double dt = Delta/fabs(val(u.y,0,0,0));




      dt *= val(cm,0,0,0);

      if (dt < dtmax) dtmax = dt;
    }}end_is_face_y()
#line 6
is_face_z(){
    if (val(u.z,0,0,0) != 0.) {
      double dt = Delta/fabs(val(u.z,0,0,0));




      dt *= val(cm,0,0,0);

      if (dt < dtmax) dtmax = dt;
    }}end_is_face_z()}end_foreach_face_generic();mpi_all_reduce_array(&dtmax,double,MPI_MIN,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 16
}else {double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);
  
#line 6
foreach_face_stencil(1,{(NonLocal[]){{"dtmax","double",(void *)&dtmax,NULL,0,'m'},{"_const_cm","double",(void *)&_const_cm,NULL,0},{"u","vector",(void *)&u,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n// #line 6 \"/home/lisergey/basilisk/src/timestep.h\"\n{is_face_x(){\n    if (val(u.x,0,0,0) != 0.) {\n      real dt = Delta/fabs(val(u.x,0,0,0));\n\n\n\n\n      dt *= _const_cm;\n\n      if (dt < dtmax) dtmax = dt;\n    }}end_is_face_x()\n// #line 6\nis_face_y(){\n    if (val(u.y,0,0,0) != 0.) {\n      real dt = Delta/fabs(val(u.y,0,0,0));\n\n\n\n\n      dt *= _const_cm;\n\n      if (dt < dtmax) dtmax = dt;\n    }}end_is_face_y()\n// #line 6\nis_face_z(){\n    if (val(u.z,0,0,0) != 0.) {\n      real dt = Delta/fabs(val(u.z,0,0,0));\n\n\n\n\n      dt *= _const_cm;\n\n      if (dt < dtmax) dtmax = dt;\n    }}end_is_face_z()}"}){_stencil_is_face_x(){
    {_stencil_val(u.x,0,0,0); {   
      _stencil_val(u.x,0,0,0);




;   




       

         
    
#line 16
}   }}end__stencil_is_face_x()
#line 6
_stencil_is_face_y(){
    {_stencil_val(u.y,0,0,0); {   
      _stencil_val(u.y,0,0,0);




;   




       

         
    
#line 16
}   }}end__stencil_is_face_y()
#line 6
_stencil_is_face_z(){
    {_stencil_val(u.z,0,0,0); {   
      _stencil_val(u.z,0,0,0);




;   




       

         
    
#line 16
}   }}end__stencil_is_face_z()}end_foreach_face_stencil();
  
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel reduction(min:dtmax)){
#line 6
foreach_face_generic(){is_face_x(){
    if (val(u.x,0,0,0) != 0.) {
      double dt = Delta/fabs(val(u.x,0,0,0));




      dt *= _const_cm;

      if (dt < dtmax) dtmax = dt;
    }}end_is_face_x()
#line 6
is_face_y(){
    if (val(u.y,0,0,0) != 0.) {
      double dt = Delta/fabs(val(u.y,0,0,0));




      dt *= _const_cm;

      if (dt < dtmax) dtmax = dt;
    }}end_is_face_y()
#line 6
is_face_z(){
    if (val(u.z,0,0,0) != 0.) {
      double dt = Delta/fabs(val(u.z,0,0,0));




      dt *= _const_cm;

      if (dt < dtmax) dtmax = dt;
    }}end_is_face_z()}end_foreach_face_generic();mpi_all_reduce_array(&dtmax,double,MPI_MIN,1);
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
#line 29 "/home/lisergey/basilisk/src/navier-stokes/centered.h"
#line 1 "./bcg.h"
#line 1 "/home/lisergey/basilisk/src/bcg.h"
#line 11 "/home/lisergey/basilisk/src/bcg.h"
void tracer_fluxes (scalar f,
      vector uf,
      vector flux,
      double dt,
              scalar src)
{





  vector  g=new_vector("g");
  gradients (((scalar[]){f,{-1}}),((vector[]) {g,{{-1},{-1},{-1}}}));




  if(!is_constant(fm.x) && !is_constant(src)){




  
#line 28
foreach_face_stencil(1,{(NonLocal[]){{"flux","vector",(void *)&flux,NULL,0},{"g","vector",(void *)&g,NULL,0},{"src","scalar",(void *)&src,NULL,0},{"f","scalar",(void *)&f,NULL,0},{"fm","vector",(void *)&fm,NULL,0},{"uf","vector",(void *)&uf,NULL,0},{"dt","double",(void *)&dt,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n// #line 28 \"/home/lisergey/basilisk/src/bcg.h\"\n{is_face_x(){ {\n\n\n\n\n\n\n\n    real un = dt*val(uf.x,0,0,0)/(val(fm.x,0,0,0)*Delta + 0.), s = sign(un);\n    int i = -(s + 1.)/2.;\n    real f2 = val(f,i,0,0) + (val(src,0,0,0) + val(src,-1,0,0))*dt/4. + s*(1. - s*un)*val(g.x,i,0,0)*Delta/2.;\n\n\n\n\n\n    if (val(fm.y,i,0,0) && val(fm.y,i,1,0)) {\n      real vn = (val(uf.y,i,0,0) + val(uf.y,i,1,0))/(val(fm.y,i,0,0) + val(fm.y,i,1,0));\n      real fyy = vn < 0. ? val(f,i,1,0) - val(f,i,0,0) : val(f,i,0,0) - val(f,i,-1,0);\n      f2 -= dt*vn*fyy/(2.*Delta);\n    }\n\n\n    if (val(fm.z,i,0,0) && val(fm.z,i,0,1)) {\n      real wn = (val(uf.z,i,0,0) + val(uf.z,i,0,1))/(val(fm.z,i,0,0) + val(fm.z,i,0,1));\n      real fzz = wn < 0. ? val(f,i,0,1) - val(f,i,0,0) : val(f,i,0,0) - val(f,i,0,-1);\n      f2 -= dt*wn*fzz/(2.*Delta);\n    }\n\n\n    val_out_(flux.x,0,0,0) = f2*val(uf.x,0,0,0);\n  }}end_is_face_x()\n// #line 28\nis_face_y(){ {\n\n\n\n\n\n\n\n    real un = dt*val(uf.y,0,0,0)/(val(fm.y,0,0,0)*Delta + 0.), s = sign(un);\n    int i = -(s + 1.)/2.;\n    real f2 = val(f,0,i,0) + (val(src,0,0,0) + val(src,0,-1,0))*dt/4. + s*(1. - s*un)*val(g.y,0,i,0)*Delta/2.;\n\n\n\n\n\n    if (val(fm.z,0,i,0) && val(fm.z,0,i,1)) {\n      real vn = (val(uf.z,0,i,0) + val(uf.z,0,i,1))/(val(fm.z,0,i,0) + val(fm.z,0,i,1));\n      real fyy = vn < 0. ? val(f,0,i,1) - val(f,0,i,0) : val(f,0,i,0) - val(f,0,i,-1);\n      f2 -= dt*vn*fyy/(2.*Delta);\n    }\n\n\n    if (val(fm.x,0,i,0) && val(fm.x,1,i,0)) {\n      real wn = (val(uf.x,0,i,0) + val(uf.x,1,i,0))/(val(fm.x,0,i,0) + val(fm.x,1,i,0));\n      real fzz = wn < 0. ? val(f,1,i,0) - val(f,0,i,0) : val(f,0,i,0) - val(f,-1,i,0);\n      f2 -= dt*wn*fzz/(2.*Delta);\n    }\n\n\n    val_out_(flux.y,0,0,0) = f2*val(uf.y,0,0,0);\n  }}end_is_face_y()\n// #line 28\nis_face_z(){ {\n\n\n\n\n\n\n\n    real un = dt*val(uf.z,0,0,0)/(val(fm.z,0,0,0)*Delta + 0.), s = sign(un);\n    int i = -(s + 1.)/2.;\n    real f2 = val(f,0,0,i) + (val(src,0,0,0) + val(src,0,0,-1))*dt/4. + s*(1. - s*un)*val(g.z,0,0,i)*Delta/2.;\n\n\n\n\n\n    if (val(fm.x,0,0,i) && val(fm.x,1,0,i)) {\n      real vn = (val(uf.x,0,0,i) + val(uf.x,1,0,i))/(val(fm.x,0,0,i) + val(fm.x,1,0,i));\n      real fyy = vn < 0. ? val(f,1,0,i) - val(f,0,0,i) : val(f,0,0,i) - val(f,-1,0,i);\n      f2 -= dt*vn*fyy/(2.*Delta);\n    }\n\n\n    if (val(fm.y,0,0,i) && val(fm.y,0,1,i)) {\n      real wn = (val(uf.y,0,0,i) + val(uf.y,0,1,i))/(val(fm.y,0,0,i) + val(fm.y,0,1,i));\n      real fzz = wn < 0. ? val(f,0,1,i) - val(f,0,0,i) : val(f,0,0,i) - val(f,0,-1,i);\n      f2 -= dt*wn*fzz/(2.*Delta);\n    }\n\n\n    val_out_(flux.z,0,0,0) = f2*val(uf.z,0,0,0);\n  }}end_is_face_z()}"}){_stencil_is_face_x(){ {        







    _stencil_val(fm.x,0,0,0);_stencil_val(uf.x,0,0,0);              
    
    _stencil_val(g.x,o_stencil,0,0); _stencil_val(src,-1,0,0);_stencil_val(src,0,0,0);_stencil_val(f, o_stencil,0,0);





_stencil_val(fm.y,o_stencil,0,0); _stencil_val(fm.y,o_stencil,1,0); {     
       _stencil_val(fm.y,o_stencil,1,0);_stencil_val(fm.y,o_stencil,0,0); _stencil_val(uf.y,o_stencil,1,0);_stencil_val(uf.y,o_stencil,0,0);         
       _stencil_val(f,o_stencil,-1,0);_stencil_val(f, o_stencil,0,0);_stencil_val(f, o_stencil,0,0); _stencil_val(f,o_stencil,1,0);
        
    }


_stencil_val(fm.z,o_stencil,0,0); _stencil_val(fm.z,o_stencil,0,1); {     
       _stencil_val(fm.z,o_stencil,0,1);_stencil_val(fm.z,o_stencil,0,0); _stencil_val(uf.z,o_stencil,0,1);_stencil_val(uf.z,o_stencil,0,0);         
       _stencil_val(f,o_stencil,0,-1);_stencil_val(f, o_stencil,0,0);_stencil_val(f, o_stencil,0,0); _stencil_val(f,o_stencil,0,1);
        
    }


_stencil_val(uf.x,0,0,0);





      


      


    
#line 58
_stencil_val_a(flux.x,0,0,0);  
  }}end__stencil_is_face_x()
#line 28
_stencil_is_face_y(){ {        







    _stencil_val(fm.y,0,0,0);_stencil_val(uf.y,0,0,0);              
    
    _stencil_val(g.y,0,o_stencil,0); _stencil_val(src,0,-1,0);_stencil_val(src,0,0,0);_stencil_val(f,0, o_stencil,0);





_stencil_val(fm.z,0,o_stencil,0); _stencil_val(fm.z,0,o_stencil,1); {     
       _stencil_val(fm.z,0,o_stencil,1);_stencil_val(fm.z,0,o_stencil,0); _stencil_val(uf.z,0,o_stencil,1);_stencil_val(uf.z,0,o_stencil,0);         
       _stencil_val(f,0,o_stencil,-1);_stencil_val(f,0, o_stencil,0);_stencil_val(f,0, o_stencil,0); _stencil_val(f,0,o_stencil,1);
        
    }


_stencil_val(fm.x,0,o_stencil,0); _stencil_val(fm.x,1,o_stencil,0); {     
       _stencil_val(fm.x,1,o_stencil,0);_stencil_val(fm.x,0,o_stencil,0); _stencil_val(uf.x,1,o_stencil,0);_stencil_val(uf.x,0,o_stencil,0);         
       _stencil_val(f,-1,o_stencil,0);_stencil_val(f,0, o_stencil,0);_stencil_val(f,0, o_stencil,0); _stencil_val(f,1,o_stencil,0);
        
    }


_stencil_val(uf.y,0,0,0);





      


      


    
#line 58
_stencil_val_a(flux.y,0,0,0);  
  }}end__stencil_is_face_y()
#line 28
_stencil_is_face_z(){ {        







    _stencil_val(fm.z,0,0,0);_stencil_val(uf.z,0,0,0);              
    
    _stencil_val(g.z,0,0,o_stencil); _stencil_val(src,0,0,-1);_stencil_val(src,0,0,0);_stencil_val(f,0,0, o_stencil);





_stencil_val(fm.x,0,0,o_stencil); _stencil_val(fm.x,1,0,o_stencil); {     
       _stencil_val(fm.x,1,0,o_stencil);_stencil_val(fm.x,0,0,o_stencil); _stencil_val(uf.x,1,0,o_stencil);_stencil_val(uf.x,0,0,o_stencil);         
       _stencil_val(f,-1,0,o_stencil);_stencil_val(f,0,0, o_stencil);_stencil_val(f,0,0, o_stencil); _stencil_val(f,1,0,o_stencil);
        
    }


_stencil_val(fm.y,0,0,o_stencil); _stencil_val(fm.y,0,1,o_stencil); {     
       _stencil_val(fm.y,0,1,o_stencil);_stencil_val(fm.y,0,0,o_stencil); _stencil_val(uf.y,0,1,o_stencil);_stencil_val(uf.y,0,0,o_stencil);         
       _stencil_val(f,0,-1,o_stencil);_stencil_val(f,0,0, o_stencil);_stencil_val(f,0,0, o_stencil); _stencil_val(f,0,1,o_stencil);
        
    }


_stencil_val(uf.z,0,0,0);





      


      


    
#line 58
_stencil_val_a(flux.z,0,0,0);  
  }}end__stencil_is_face_z()}end_foreach_face_stencil();{
#line 28
foreach_face_generic(){is_face_x(){ {







    double un = dt*val(uf.x,0,0,0)/(val(fm.x,0,0,0)*Delta + 0.), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = val(f,i,0,0) + (val(src,0,0,0) + val(src,-1,0,0))*dt/4. + s*(1. - s*un)*val(g.x,i,0,0)*Delta/2.;





    if (val(fm.y,i,0,0) && val(fm.y,i,1,0)) {
      double vn = (val(uf.y,i,0,0) + val(uf.y,i,1,0))/(val(fm.y,i,0,0) + val(fm.y,i,1,0));
      double fyy = vn < 0. ? val(f,i,1,0) - val(f,i,0,0) : val(f,i,0,0) - val(f,i,-1,0);
      f2 -= dt*vn*fyy/(2.*Delta);
    }


    if (val(fm.z,i,0,0) && val(fm.z,i,0,1)) {
      double wn = (val(uf.z,i,0,0) + val(uf.z,i,0,1))/(val(fm.z,i,0,0) + val(fm.z,i,0,1));
      double fzz = wn < 0. ? val(f,i,0,1) - val(f,i,0,0) : val(f,i,0,0) - val(f,i,0,-1);
      f2 -= dt*wn*fzz/(2.*Delta);
    }


    val(flux.x,0,0,0) = f2*val(uf.x,0,0,0);
  }}end_is_face_x()
#line 28
is_face_y(){ {







    double un = dt*val(uf.y,0,0,0)/(val(fm.y,0,0,0)*Delta + 0.), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = val(f,0,i,0) + (val(src,0,0,0) + val(src,0,-1,0))*dt/4. + s*(1. - s*un)*val(g.y,0,i,0)*Delta/2.;





    if (val(fm.z,0,i,0) && val(fm.z,0,i,1)) {
      double vn = (val(uf.z,0,i,0) + val(uf.z,0,i,1))/(val(fm.z,0,i,0) + val(fm.z,0,i,1));
      double fyy = vn < 0. ? val(f,0,i,1) - val(f,0,i,0) : val(f,0,i,0) - val(f,0,i,-1);
      f2 -= dt*vn*fyy/(2.*Delta);
    }


    if (val(fm.x,0,i,0) && val(fm.x,1,i,0)) {
      double wn = (val(uf.x,0,i,0) + val(uf.x,1,i,0))/(val(fm.x,0,i,0) + val(fm.x,1,i,0));
      double fzz = wn < 0. ? val(f,1,i,0) - val(f,0,i,0) : val(f,0,i,0) - val(f,-1,i,0);
      f2 -= dt*wn*fzz/(2.*Delta);
    }


    val(flux.y,0,0,0) = f2*val(uf.y,0,0,0);
  }}end_is_face_y()
#line 28
is_face_z(){ {







    double un = dt*val(uf.z,0,0,0)/(val(fm.z,0,0,0)*Delta + 0.), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = val(f,0,0,i) + (val(src,0,0,0) + val(src,0,0,-1))*dt/4. + s*(1. - s*un)*val(g.z,0,0,i)*Delta/2.;





    if (val(fm.x,0,0,i) && val(fm.x,1,0,i)) {
      double vn = (val(uf.x,0,0,i) + val(uf.x,1,0,i))/(val(fm.x,0,0,i) + val(fm.x,1,0,i));
      double fyy = vn < 0. ? val(f,1,0,i) - val(f,0,0,i) : val(f,0,0,i) - val(f,-1,0,i);
      f2 -= dt*vn*fyy/(2.*Delta);
    }


    if (val(fm.y,0,0,i) && val(fm.y,0,1,i)) {
      double wn = (val(uf.y,0,0,i) + val(uf.y,0,1,i))/(val(fm.y,0,0,i) + val(fm.y,0,1,i));
      double fzz = wn < 0. ? val(f,0,1,i) - val(f,0,0,i) : val(f,0,0,i) - val(f,0,-1,i);
      f2 -= dt*wn*fzz/(2.*Delta);
    }


    val(flux.z,0,0,0) = f2*val(uf.z,0,0,0);
  }}end_is_face_z()}end_foreach_face_generic();}}else if(is_constant(fm.x) && !is_constant(src)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);




  
#line 28
foreach_face_stencil(1,{(NonLocal[]){{"flux","vector",(void *)&flux,NULL,0},{"g","vector",(void *)&g,NULL,0},{"src","scalar",(void *)&src,NULL,0},{"f","scalar",(void *)&f,NULL,0},{"_const_fm","_coord",(void *)&_const_fm,NULL,0},{"uf","vector",(void *)&uf,NULL,0},{"dt","double",(void *)&dt,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n// #line 28 \"/home/lisergey/basilisk/src/bcg.h\"\n{is_face_x(){ {\n\n\n\n\n\n\n\n    real un = dt*val(uf.x,0,0,0)/(_const_fm.x*Delta + 0.), s = sign(un);\n    int i = -(s + 1.)/2.;\n    real f2 = val(f,i,0,0) + (val(src,0,0,0) + val(src,-1,0,0))*dt/4. + s*(1. - s*un)*val(g.x,i,0,0)*Delta/2.;\n\n\n\n\n\n    if (_const_fm.y && _const_fm.y) {\n      real vn = (val(uf.y,i,0,0) + val(uf.y,i,1,0))/(_const_fm.y + _const_fm.y);\n      real fyy = vn < 0. ? val(f,i,1,0) - val(f,i,0,0) : val(f,i,0,0) - val(f,i,-1,0);\n      f2 -= dt*vn*fyy/(2.*Delta);\n    }\n\n\n    if (_const_fm.z && _const_fm.z) {\n      real wn = (val(uf.z,i,0,0) + val(uf.z,i,0,1))/(_const_fm.z + _const_fm.z);\n      real fzz = wn < 0. ? val(f,i,0,1) - val(f,i,0,0) : val(f,i,0,0) - val(f,i,0,-1);\n      f2 -= dt*wn*fzz/(2.*Delta);\n    }\n\n\n    val_out_(flux.x,0,0,0) = f2*val(uf.x,0,0,0);\n  }}end_is_face_x()\n// #line 28\nis_face_y(){ {\n\n\n\n\n\n\n\n    real un = dt*val(uf.y,0,0,0)/(_const_fm.y*Delta + 0.), s = sign(un);\n    int i = -(s + 1.)/2.;\n    real f2 = val(f,0,i,0) + (val(src,0,0,0) + val(src,0,-1,0))*dt/4. + s*(1. - s*un)*val(g.y,0,i,0)*Delta/2.;\n\n\n\n\n\n    if (_const_fm.z && _const_fm.z) {\n      real vn = (val(uf.z,0,i,0) + val(uf.z,0,i,1))/(_const_fm.z + _const_fm.z);\n      real fyy = vn < 0. ? val(f,0,i,1) - val(f,0,i,0) : val(f,0,i,0) - val(f,0,i,-1);\n      f2 -= dt*vn*fyy/(2.*Delta);\n    }\n\n\n    if (_const_fm.x && _const_fm.x) {\n      real wn = (val(uf.x,0,i,0) + val(uf.x,1,i,0))/(_const_fm.x + _const_fm.x);\n      real fzz = wn < 0. ? val(f,1,i,0) - val(f,0,i,0) : val(f,0,i,0) - val(f,-1,i,0);\n      f2 -= dt*wn*fzz/(2.*Delta);\n    }\n\n\n    val_out_(flux.y,0,0,0) = f2*val(uf.y,0,0,0);\n  }}end_is_face_y()\n// #line 28\nis_face_z(){ {\n\n\n\n\n\n\n\n    real un = dt*val(uf.z,0,0,0)/(_const_fm.z*Delta + 0.), s = sign(un);\n    int i = -(s + 1.)/2.;\n    real f2 = val(f,0,0,i) + (val(src,0,0,0) + val(src,0,0,-1))*dt/4. + s*(1. - s*un)*val(g.z,0,0,i)*Delta/2.;\n\n\n\n\n\n    if (_const_fm.x && _const_fm.x) {\n      real vn = (val(uf.x,0,0,i) + val(uf.x,1,0,i))/(_const_fm.x + _const_fm.x);\n      real fyy = vn < 0. ? val(f,1,0,i) - val(f,0,0,i) : val(f,0,0,i) - val(f,-1,0,i);\n      f2 -= dt*vn*fyy/(2.*Delta);\n    }\n\n\n    if (_const_fm.y && _const_fm.y) {\n      real wn = (val(uf.y,0,0,i) + val(uf.y,0,1,i))/(_const_fm.y + _const_fm.y);\n      real fzz = wn < 0. ? val(f,0,1,i) - val(f,0,0,i) : val(f,0,0,i) - val(f,0,-1,i);\n      f2 -= dt*wn*fzz/(2.*Delta);\n    }\n\n\n    val_out_(flux.z,0,0,0) = f2*val(uf.z,0,0,0);\n  }}end_is_face_z()}"}){_stencil_is_face_x(){ {







;_stencil_val(uf.x,0,0,0);              
    
    _stencil_val(g.x,o_stencil,0,0); _stencil_val(src,-1,0,0);_stencil_val(src,0,0,0);_stencil_val(f, o_stencil,0,0);





;; {
;; _stencil_val(uf.y,o_stencil,1,0);_stencil_val(uf.y,o_stencil,0,0);         
       _stencil_val(f,o_stencil,-1,0);_stencil_val(f, o_stencil,0,0);_stencil_val(f, o_stencil,0,0); _stencil_val(f,o_stencil,1,0);
        
    }


;; {
;; _stencil_val(uf.z,o_stencil,0,1);_stencil_val(uf.z,o_stencil,0,0);         
       _stencil_val(f,o_stencil,0,-1);_stencil_val(f, o_stencil,0,0);_stencil_val(f, o_stencil,0,0); _stencil_val(f,o_stencil,0,1);
        
    }


_stencil_val(uf.x,0,0,0);





      


      


    
#line 58
_stencil_val_a(flux.x,0,0,0);  
  }}end__stencil_is_face_x()
#line 28
_stencil_is_face_y(){ {







;_stencil_val(uf.y,0,0,0);              
    
    _stencil_val(g.y,0,o_stencil,0); _stencil_val(src,0,-1,0);_stencil_val(src,0,0,0);_stencil_val(f,0, o_stencil,0);





;; {
;; _stencil_val(uf.z,0,o_stencil,1);_stencil_val(uf.z,0,o_stencil,0);         
       _stencil_val(f,0,o_stencil,-1);_stencil_val(f,0, o_stencil,0);_stencil_val(f,0, o_stencil,0); _stencil_val(f,0,o_stencil,1);
        
    }


;; {
;; _stencil_val(uf.x,1,o_stencil,0);_stencil_val(uf.x,0,o_stencil,0);         
       _stencil_val(f,-1,o_stencil,0);_stencil_val(f,0, o_stencil,0);_stencil_val(f,0, o_stencil,0); _stencil_val(f,1,o_stencil,0);
        
    }


_stencil_val(uf.y,0,0,0);





      


      


    
#line 58
_stencil_val_a(flux.y,0,0,0);  
  }}end__stencil_is_face_y()
#line 28
_stencil_is_face_z(){ {







;_stencil_val(uf.z,0,0,0);              
    
    _stencil_val(g.z,0,0,o_stencil); _stencil_val(src,0,0,-1);_stencil_val(src,0,0,0);_stencil_val(f,0,0, o_stencil);





;; {
;; _stencil_val(uf.x,1,0,o_stencil);_stencil_val(uf.x,0,0,o_stencil);         
       _stencil_val(f,-1,0,o_stencil);_stencil_val(f,0,0, o_stencil);_stencil_val(f,0,0, o_stencil); _stencil_val(f,1,0,o_stencil);
        
    }


;; {
;; _stencil_val(uf.y,0,1,o_stencil);_stencil_val(uf.y,0,0,o_stencil);         
       _stencil_val(f,0,-1,o_stencil);_stencil_val(f,0,0, o_stencil);_stencil_val(f,0,0, o_stencil); _stencil_val(f,0,1,o_stencil);
        
    }


_stencil_val(uf.z,0,0,0);





      


      


    
#line 58
_stencil_val_a(flux.z,0,0,0);  
  }}end__stencil_is_face_z()}end_foreach_face_stencil();




  {
#line 28
foreach_face_generic(){is_face_x(){ {







    double un = dt*val(uf.x,0,0,0)/(_const_fm.x*Delta + 0.), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = val(f,i,0,0) + (val(src,0,0,0) + val(src,-1,0,0))*dt/4. + s*(1. - s*un)*val(g.x,i,0,0)*Delta/2.;





    if (_const_fm.y && _const_fm.y) {
      double vn = (val(uf.y,i,0,0) + val(uf.y,i,1,0))/(_const_fm.y + _const_fm.y);
      double fyy = vn < 0. ? val(f,i,1,0) - val(f,i,0,0) : val(f,i,0,0) - val(f,i,-1,0);
      f2 -= dt*vn*fyy/(2.*Delta);
    }


    if (_const_fm.z && _const_fm.z) {
      double wn = (val(uf.z,i,0,0) + val(uf.z,i,0,1))/(_const_fm.z + _const_fm.z);
      double fzz = wn < 0. ? val(f,i,0,1) - val(f,i,0,0) : val(f,i,0,0) - val(f,i,0,-1);
      f2 -= dt*wn*fzz/(2.*Delta);
    }


    val(flux.x,0,0,0) = f2*val(uf.x,0,0,0);
  }}end_is_face_x()
#line 28
is_face_y(){ {







    double un = dt*val(uf.y,0,0,0)/(_const_fm.y*Delta + 0.), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = val(f,0,i,0) + (val(src,0,0,0) + val(src,0,-1,0))*dt/4. + s*(1. - s*un)*val(g.y,0,i,0)*Delta/2.;





    if (_const_fm.z && _const_fm.z) {
      double vn = (val(uf.z,0,i,0) + val(uf.z,0,i,1))/(_const_fm.z + _const_fm.z);
      double fyy = vn < 0. ? val(f,0,i,1) - val(f,0,i,0) : val(f,0,i,0) - val(f,0,i,-1);
      f2 -= dt*vn*fyy/(2.*Delta);
    }


    if (_const_fm.x && _const_fm.x) {
      double wn = (val(uf.x,0,i,0) + val(uf.x,1,i,0))/(_const_fm.x + _const_fm.x);
      double fzz = wn < 0. ? val(f,1,i,0) - val(f,0,i,0) : val(f,0,i,0) - val(f,-1,i,0);
      f2 -= dt*wn*fzz/(2.*Delta);
    }


    val(flux.y,0,0,0) = f2*val(uf.y,0,0,0);
  }}end_is_face_y()
#line 28
is_face_z(){ {







    double un = dt*val(uf.z,0,0,0)/(_const_fm.z*Delta + 0.), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = val(f,0,0,i) + (val(src,0,0,0) + val(src,0,0,-1))*dt/4. + s*(1. - s*un)*val(g.z,0,0,i)*Delta/2.;





    if (_const_fm.x && _const_fm.x) {
      double vn = (val(uf.x,0,0,i) + val(uf.x,1,0,i))/(_const_fm.x + _const_fm.x);
      double fyy = vn < 0. ? val(f,1,0,i) - val(f,0,0,i) : val(f,0,0,i) - val(f,-1,0,i);
      f2 -= dt*vn*fyy/(2.*Delta);
    }


    if (_const_fm.y && _const_fm.y) {
      double wn = (val(uf.y,0,0,i) + val(uf.y,0,1,i))/(_const_fm.y + _const_fm.y);
      double fzz = wn < 0. ? val(f,0,1,i) - val(f,0,0,i) : val(f,0,0,i) - val(f,0,-1,i);
      f2 -= dt*wn*fzz/(2.*Delta);
    }


    val(flux.z,0,0,0) = f2*val(uf.z,0,0,0);
  }}end_is_face_z()}end_foreach_face_generic();}}else if(!is_constant(fm.x) && is_constant(src)){double _const_src=_constant[src.i-_NVARMAX];NOT_UNUSED(_const_src);




  
#line 28
foreach_face_stencil(1,{(NonLocal[]){{"flux","vector",(void *)&flux,NULL,0},{"g","vector",(void *)&g,NULL,0},{"_const_src","double",(void *)&_const_src,NULL,0},{"f","scalar",(void *)&f,NULL,0},{"fm","vector",(void *)&fm,NULL,0},{"uf","vector",(void *)&uf,NULL,0},{"dt","double",(void *)&dt,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n// #line 28 \"/home/lisergey/basilisk/src/bcg.h\"\n{is_face_x(){ {\n\n\n\n\n\n\n\n    real un = dt*val(uf.x,0,0,0)/(val(fm.x,0,0,0)*Delta + 0.), s = sign(un);\n    int i = -(s + 1.)/2.;\n    real f2 = val(f,i,0,0) + (_const_src + _const_src)*dt/4. + s*(1. - s*un)*val(g.x,i,0,0)*Delta/2.;\n\n\n\n\n\n    if (val(fm.y,i,0,0) && val(fm.y,i,1,0)) {\n      real vn = (val(uf.y,i,0,0) + val(uf.y,i,1,0))/(val(fm.y,i,0,0) + val(fm.y,i,1,0));\n      real fyy = vn < 0. ? val(f,i,1,0) - val(f,i,0,0) : val(f,i,0,0) - val(f,i,-1,0);\n      f2 -= dt*vn*fyy/(2.*Delta);\n    }\n\n\n    if (val(fm.z,i,0,0) && val(fm.z,i,0,1)) {\n      real wn = (val(uf.z,i,0,0) + val(uf.z,i,0,1))/(val(fm.z,i,0,0) + val(fm.z,i,0,1));\n      real fzz = wn < 0. ? val(f,i,0,1) - val(f,i,0,0) : val(f,i,0,0) - val(f,i,0,-1);\n      f2 -= dt*wn*fzz/(2.*Delta);\n    }\n\n\n    val_out_(flux.x,0,0,0) = f2*val(uf.x,0,0,0);\n  }}end_is_face_x()\n// #line 28\nis_face_y(){ {\n\n\n\n\n\n\n\n    real un = dt*val(uf.y,0,0,0)/(val(fm.y,0,0,0)*Delta + 0.), s = sign(un);\n    int i = -(s + 1.)/2.;\n    real f2 = val(f,0,i,0) + (_const_src + _const_src)*dt/4. + s*(1. - s*un)*val(g.y,0,i,0)*Delta/2.;\n\n\n\n\n\n    if (val(fm.z,0,i,0) && val(fm.z,0,i,1)) {\n      real vn = (val(uf.z,0,i,0) + val(uf.z,0,i,1))/(val(fm.z,0,i,0) + val(fm.z,0,i,1));\n      real fyy = vn < 0. ? val(f,0,i,1) - val(f,0,i,0) : val(f,0,i,0) - val(f,0,i,-1);\n      f2 -= dt*vn*fyy/(2.*Delta);\n    }\n\n\n    if (val(fm.x,0,i,0) && val(fm.x,1,i,0)) {\n      real wn = (val(uf.x,0,i,0) + val(uf.x,1,i,0))/(val(fm.x,0,i,0) + val(fm.x,1,i,0));\n      real fzz = wn < 0. ? val(f,1,i,0) - val(f,0,i,0) : val(f,0,i,0) - val(f,-1,i,0);\n      f2 -= dt*wn*fzz/(2.*Delta);\n    }\n\n\n    val_out_(flux.y,0,0,0) = f2*val(uf.y,0,0,0);\n  }}end_is_face_y()\n// #line 28\nis_face_z(){ {\n\n\n\n\n\n\n\n    real un = dt*val(uf.z,0,0,0)/(val(fm.z,0,0,0)*Delta + 0.), s = sign(un);\n    int i = -(s + 1.)/2.;\n    real f2 = val(f,0,0,i) + (_const_src + _const_src)*dt/4. + s*(1. - s*un)*val(g.z,0,0,i)*Delta/2.;\n\n\n\n\n\n    if (val(fm.x,0,0,i) && val(fm.x,1,0,i)) {\n      real vn = (val(uf.x,0,0,i) + val(uf.x,1,0,i))/(val(fm.x,0,0,i) + val(fm.x,1,0,i));\n      real fyy = vn < 0. ? val(f,1,0,i) - val(f,0,0,i) : val(f,0,0,i) - val(f,-1,0,i);\n      f2 -= dt*vn*fyy/(2.*Delta);\n    }\n\n\n    if (val(fm.y,0,0,i) && val(fm.y,0,1,i)) {\n      real wn = (val(uf.y,0,0,i) + val(uf.y,0,1,i))/(val(fm.y,0,0,i) + val(fm.y,0,1,i));\n      real fzz = wn < 0. ? val(f,0,1,i) - val(f,0,0,i) : val(f,0,0,i) - val(f,0,-1,i);\n      f2 -= dt*wn*fzz/(2.*Delta);\n    }\n\n\n    val_out_(flux.z,0,0,0) = f2*val(uf.z,0,0,0);\n  }}end_is_face_z()}"}){_stencil_is_face_x(){ {        







    _stencil_val(fm.x,0,0,0);_stencil_val(uf.x,0,0,0);              
    
    _stencil_val(g.x,o_stencil,0,0);;;_stencil_val(f, o_stencil,0,0);





_stencil_val(fm.y,o_stencil,0,0); _stencil_val(fm.y,o_stencil,1,0); {     
       _stencil_val(fm.y,o_stencil,1,0);_stencil_val(fm.y,o_stencil,0,0); _stencil_val(uf.y,o_stencil,1,0);_stencil_val(uf.y,o_stencil,0,0);         
       _stencil_val(f,o_stencil,-1,0);_stencil_val(f, o_stencil,0,0);_stencil_val(f, o_stencil,0,0); _stencil_val(f,o_stencil,1,0);
        
    }


_stencil_val(fm.z,o_stencil,0,0); _stencil_val(fm.z,o_stencil,0,1); {     
       _stencil_val(fm.z,o_stencil,0,1);_stencil_val(fm.z,o_stencil,0,0); _stencil_val(uf.z,o_stencil,0,1);_stencil_val(uf.z,o_stencil,0,0);         
       _stencil_val(f,o_stencil,0,-1);_stencil_val(f, o_stencil,0,0);_stencil_val(f, o_stencil,0,0); _stencil_val(f,o_stencil,0,1);
        
    }


_stencil_val(uf.x,0,0,0);





      


      


    
#line 58
_stencil_val_a(flux.x,0,0,0);  
  }}end__stencil_is_face_x()
#line 28
_stencil_is_face_y(){ {        







    _stencil_val(fm.y,0,0,0);_stencil_val(uf.y,0,0,0);              
    
    _stencil_val(g.y,0,o_stencil,0);;;_stencil_val(f,0, o_stencil,0);





_stencil_val(fm.z,0,o_stencil,0); _stencil_val(fm.z,0,o_stencil,1); {     
       _stencil_val(fm.z,0,o_stencil,1);_stencil_val(fm.z,0,o_stencil,0); _stencil_val(uf.z,0,o_stencil,1);_stencil_val(uf.z,0,o_stencil,0);         
       _stencil_val(f,0,o_stencil,-1);_stencil_val(f,0, o_stencil,0);_stencil_val(f,0, o_stencil,0); _stencil_val(f,0,o_stencil,1);
        
    }


_stencil_val(fm.x,0,o_stencil,0); _stencil_val(fm.x,1,o_stencil,0); {     
       _stencil_val(fm.x,1,o_stencil,0);_stencil_val(fm.x,0,o_stencil,0); _stencil_val(uf.x,1,o_stencil,0);_stencil_val(uf.x,0,o_stencil,0);         
       _stencil_val(f,-1,o_stencil,0);_stencil_val(f,0, o_stencil,0);_stencil_val(f,0, o_stencil,0); _stencil_val(f,1,o_stencil,0);
        
    }


_stencil_val(uf.y,0,0,0);





      


      


    
#line 58
_stencil_val_a(flux.y,0,0,0);  
  }}end__stencil_is_face_y()
#line 28
_stencil_is_face_z(){ {        







    _stencil_val(fm.z,0,0,0);_stencil_val(uf.z,0,0,0);              
    
    _stencil_val(g.z,0,0,o_stencil);;;_stencil_val(f,0,0, o_stencil);





_stencil_val(fm.x,0,0,o_stencil); _stencil_val(fm.x,1,0,o_stencil); {     
       _stencil_val(fm.x,1,0,o_stencil);_stencil_val(fm.x,0,0,o_stencil); _stencil_val(uf.x,1,0,o_stencil);_stencil_val(uf.x,0,0,o_stencil);         
       _stencil_val(f,-1,0,o_stencil);_stencil_val(f,0,0, o_stencil);_stencil_val(f,0,0, o_stencil); _stencil_val(f,1,0,o_stencil);
        
    }


_stencil_val(fm.y,0,0,o_stencil); _stencil_val(fm.y,0,1,o_stencil); {     
       _stencil_val(fm.y,0,1,o_stencil);_stencil_val(fm.y,0,0,o_stencil); _stencil_val(uf.y,0,1,o_stencil);_stencil_val(uf.y,0,0,o_stencil);         
       _stencil_val(f,0,-1,o_stencil);_stencil_val(f,0,0, o_stencil);_stencil_val(f,0,0, o_stencil); _stencil_val(f,0,1,o_stencil);
        
    }


_stencil_val(uf.z,0,0,0);





      


      


    
#line 58
_stencil_val_a(flux.z,0,0,0);  
  }}end__stencil_is_face_z()}end_foreach_face_stencil();




  {
#line 28
foreach_face_generic(){is_face_x(){ {







    double un = dt*val(uf.x,0,0,0)/(val(fm.x,0,0,0)*Delta + 0.), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = val(f,i,0,0) + (_const_src + _const_src)*dt/4. + s*(1. - s*un)*val(g.x,i,0,0)*Delta/2.;





    if (val(fm.y,i,0,0) && val(fm.y,i,1,0)) {
      double vn = (val(uf.y,i,0,0) + val(uf.y,i,1,0))/(val(fm.y,i,0,0) + val(fm.y,i,1,0));
      double fyy = vn < 0. ? val(f,i,1,0) - val(f,i,0,0) : val(f,i,0,0) - val(f,i,-1,0);
      f2 -= dt*vn*fyy/(2.*Delta);
    }


    if (val(fm.z,i,0,0) && val(fm.z,i,0,1)) {
      double wn = (val(uf.z,i,0,0) + val(uf.z,i,0,1))/(val(fm.z,i,0,0) + val(fm.z,i,0,1));
      double fzz = wn < 0. ? val(f,i,0,1) - val(f,i,0,0) : val(f,i,0,0) - val(f,i,0,-1);
      f2 -= dt*wn*fzz/(2.*Delta);
    }


    val(flux.x,0,0,0) = f2*val(uf.x,0,0,0);
  }}end_is_face_x()
#line 28
is_face_y(){ {







    double un = dt*val(uf.y,0,0,0)/(val(fm.y,0,0,0)*Delta + 0.), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = val(f,0,i,0) + (_const_src + _const_src)*dt/4. + s*(1. - s*un)*val(g.y,0,i,0)*Delta/2.;





    if (val(fm.z,0,i,0) && val(fm.z,0,i,1)) {
      double vn = (val(uf.z,0,i,0) + val(uf.z,0,i,1))/(val(fm.z,0,i,0) + val(fm.z,0,i,1));
      double fyy = vn < 0. ? val(f,0,i,1) - val(f,0,i,0) : val(f,0,i,0) - val(f,0,i,-1);
      f2 -= dt*vn*fyy/(2.*Delta);
    }


    if (val(fm.x,0,i,0) && val(fm.x,1,i,0)) {
      double wn = (val(uf.x,0,i,0) + val(uf.x,1,i,0))/(val(fm.x,0,i,0) + val(fm.x,1,i,0));
      double fzz = wn < 0. ? val(f,1,i,0) - val(f,0,i,0) : val(f,0,i,0) - val(f,-1,i,0);
      f2 -= dt*wn*fzz/(2.*Delta);
    }


    val(flux.y,0,0,0) = f2*val(uf.y,0,0,0);
  }}end_is_face_y()
#line 28
is_face_z(){ {







    double un = dt*val(uf.z,0,0,0)/(val(fm.z,0,0,0)*Delta + 0.), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = val(f,0,0,i) + (_const_src + _const_src)*dt/4. + s*(1. - s*un)*val(g.z,0,0,i)*Delta/2.;





    if (val(fm.x,0,0,i) && val(fm.x,1,0,i)) {
      double vn = (val(uf.x,0,0,i) + val(uf.x,1,0,i))/(val(fm.x,0,0,i) + val(fm.x,1,0,i));
      double fyy = vn < 0. ? val(f,1,0,i) - val(f,0,0,i) : val(f,0,0,i) - val(f,-1,0,i);
      f2 -= dt*vn*fyy/(2.*Delta);
    }


    if (val(fm.y,0,0,i) && val(fm.y,0,1,i)) {
      double wn = (val(uf.y,0,0,i) + val(uf.y,0,1,i))/(val(fm.y,0,0,i) + val(fm.y,0,1,i));
      double fzz = wn < 0. ? val(f,0,1,i) - val(f,0,0,i) : val(f,0,0,i) - val(f,0,-1,i);
      f2 -= dt*wn*fzz/(2.*Delta);
    }


    val(flux.z,0,0,0) = f2*val(uf.z,0,0,0);
  }}end_is_face_z()}end_foreach_face_generic();}}else {_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_src=_constant[src.i-_NVARMAX];NOT_UNUSED(_const_src);




  
#line 28
foreach_face_stencil(1,{(NonLocal[]){{"flux","vector",(void *)&flux,NULL,0},{"g","vector",(void *)&g,NULL,0},{"_const_src","double",(void *)&_const_src,NULL,0},{"f","scalar",(void *)&f,NULL,0},{"_const_fm","_coord",(void *)&_const_fm,NULL,0},{"uf","vector",(void *)&uf,NULL,0},{"dt","double",(void *)&dt,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n// #line 28 \"/home/lisergey/basilisk/src/bcg.h\"\n{is_face_x(){ {\n\n\n\n\n\n\n\n    real un = dt*val(uf.x,0,0,0)/(_const_fm.x*Delta + 0.), s = sign(un);\n    int i = -(s + 1.)/2.;\n    real f2 = val(f,i,0,0) + (_const_src + _const_src)*dt/4. + s*(1. - s*un)*val(g.x,i,0,0)*Delta/2.;\n\n\n\n\n\n    if (_const_fm.y && _const_fm.y) {\n      real vn = (val(uf.y,i,0,0) + val(uf.y,i,1,0))/(_const_fm.y + _const_fm.y);\n      real fyy = vn < 0. ? val(f,i,1,0) - val(f,i,0,0) : val(f,i,0,0) - val(f,i,-1,0);\n      f2 -= dt*vn*fyy/(2.*Delta);\n    }\n\n\n    if (_const_fm.z && _const_fm.z) {\n      real wn = (val(uf.z,i,0,0) + val(uf.z,i,0,1))/(_const_fm.z + _const_fm.z);\n      real fzz = wn < 0. ? val(f,i,0,1) - val(f,i,0,0) : val(f,i,0,0) - val(f,i,0,-1);\n      f2 -= dt*wn*fzz/(2.*Delta);\n    }\n\n\n    val_out_(flux.x,0,0,0) = f2*val(uf.x,0,0,0);\n  }}end_is_face_x()\n// #line 28\nis_face_y(){ {\n\n\n\n\n\n\n\n    real un = dt*val(uf.y,0,0,0)/(_const_fm.y*Delta + 0.), s = sign(un);\n    int i = -(s + 1.)/2.;\n    real f2 = val(f,0,i,0) + (_const_src + _const_src)*dt/4. + s*(1. - s*un)*val(g.y,0,i,0)*Delta/2.;\n\n\n\n\n\n    if (_const_fm.z && _const_fm.z) {\n      real vn = (val(uf.z,0,i,0) + val(uf.z,0,i,1))/(_const_fm.z + _const_fm.z);\n      real fyy = vn < 0. ? val(f,0,i,1) - val(f,0,i,0) : val(f,0,i,0) - val(f,0,i,-1);\n      f2 -= dt*vn*fyy/(2.*Delta);\n    }\n\n\n    if (_const_fm.x && _const_fm.x) {\n      real wn = (val(uf.x,0,i,0) + val(uf.x,1,i,0))/(_const_fm.x + _const_fm.x);\n      real fzz = wn < 0. ? val(f,1,i,0) - val(f,0,i,0) : val(f,0,i,0) - val(f,-1,i,0);\n      f2 -= dt*wn*fzz/(2.*Delta);\n    }\n\n\n    val_out_(flux.y,0,0,0) = f2*val(uf.y,0,0,0);\n  }}end_is_face_y()\n// #line 28\nis_face_z(){ {\n\n\n\n\n\n\n\n    real un = dt*val(uf.z,0,0,0)/(_const_fm.z*Delta + 0.), s = sign(un);\n    int i = -(s + 1.)/2.;\n    real f2 = val(f,0,0,i) + (_const_src + _const_src)*dt/4. + s*(1. - s*un)*val(g.z,0,0,i)*Delta/2.;\n\n\n\n\n\n    if (_const_fm.x && _const_fm.x) {\n      real vn = (val(uf.x,0,0,i) + val(uf.x,1,0,i))/(_const_fm.x + _const_fm.x);\n      real fyy = vn < 0. ? val(f,1,0,i) - val(f,0,0,i) : val(f,0,0,i) - val(f,-1,0,i);\n      f2 -= dt*vn*fyy/(2.*Delta);\n    }\n\n\n    if (_const_fm.y && _const_fm.y) {\n      real wn = (val(uf.y,0,0,i) + val(uf.y,0,1,i))/(_const_fm.y + _const_fm.y);\n      real fzz = wn < 0. ? val(f,0,1,i) - val(f,0,0,i) : val(f,0,0,i) - val(f,0,-1,i);\n      f2 -= dt*wn*fzz/(2.*Delta);\n    }\n\n\n    val_out_(flux.z,0,0,0) = f2*val(uf.z,0,0,0);\n  }}end_is_face_z()}"}){_stencil_is_face_x(){ {







;_stencil_val(uf.x,0,0,0);              
    
    _stencil_val(g.x,o_stencil,0,0);;;_stencil_val(f, o_stencil,0,0);





;; {
;; _stencil_val(uf.y,o_stencil,1,0);_stencil_val(uf.y,o_stencil,0,0);         
       _stencil_val(f,o_stencil,-1,0);_stencil_val(f, o_stencil,0,0);_stencil_val(f, o_stencil,0,0); _stencil_val(f,o_stencil,1,0);
        
    }


;; {
;; _stencil_val(uf.z,o_stencil,0,1);_stencil_val(uf.z,o_stencil,0,0);         
       _stencil_val(f,o_stencil,0,-1);_stencil_val(f, o_stencil,0,0);_stencil_val(f, o_stencil,0,0); _stencil_val(f,o_stencil,0,1);
        
    }


_stencil_val(uf.x,0,0,0);





      


      


    
#line 58
_stencil_val_a(flux.x,0,0,0);  
  }}end__stencil_is_face_x()
#line 28
_stencil_is_face_y(){ {







;_stencil_val(uf.y,0,0,0);              
    
    _stencil_val(g.y,0,o_stencil,0);;;_stencil_val(f,0, o_stencil,0);





;; {
;; _stencil_val(uf.z,0,o_stencil,1);_stencil_val(uf.z,0,o_stencil,0);         
       _stencil_val(f,0,o_stencil,-1);_stencil_val(f,0, o_stencil,0);_stencil_val(f,0, o_stencil,0); _stencil_val(f,0,o_stencil,1);
        
    }


;; {
;; _stencil_val(uf.x,1,o_stencil,0);_stencil_val(uf.x,0,o_stencil,0);         
       _stencil_val(f,-1,o_stencil,0);_stencil_val(f,0, o_stencil,0);_stencil_val(f,0, o_stencil,0); _stencil_val(f,1,o_stencil,0);
        
    }


_stencil_val(uf.y,0,0,0);





      


      


    
#line 58
_stencil_val_a(flux.y,0,0,0);  
  }}end__stencil_is_face_y()
#line 28
_stencil_is_face_z(){ {







;_stencil_val(uf.z,0,0,0);              
    
    _stencil_val(g.z,0,0,o_stencil);;;_stencil_val(f,0,0, o_stencil);





;; {
;; _stencil_val(uf.x,1,0,o_stencil);_stencil_val(uf.x,0,0,o_stencil);         
       _stencil_val(f,-1,0,o_stencil);_stencil_val(f,0,0, o_stencil);_stencil_val(f,0,0, o_stencil); _stencil_val(f,1,0,o_stencil);
        
    }


;; {
;; _stencil_val(uf.y,0,1,o_stencil);_stencil_val(uf.y,0,0,o_stencil);         
       _stencil_val(f,0,-1,o_stencil);_stencil_val(f,0,0, o_stencil);_stencil_val(f,0,0, o_stencil); _stencil_val(f,0,1,o_stencil);
        
    }


_stencil_val(uf.z,0,0,0);





      


      


    
#line 58
_stencil_val_a(flux.z,0,0,0);  
  }}end__stencil_is_face_z()}end_foreach_face_stencil();




  {
#line 28
foreach_face_generic(){is_face_x(){ {







    double un = dt*val(uf.x,0,0,0)/(_const_fm.x*Delta + 0.), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = val(f,i,0,0) + (_const_src + _const_src)*dt/4. + s*(1. - s*un)*val(g.x,i,0,0)*Delta/2.;





    if (_const_fm.y && _const_fm.y) {
      double vn = (val(uf.y,i,0,0) + val(uf.y,i,1,0))/(_const_fm.y + _const_fm.y);
      double fyy = vn < 0. ? val(f,i,1,0) - val(f,i,0,0) : val(f,i,0,0) - val(f,i,-1,0);
      f2 -= dt*vn*fyy/(2.*Delta);
    }


    if (_const_fm.z && _const_fm.z) {
      double wn = (val(uf.z,i,0,0) + val(uf.z,i,0,1))/(_const_fm.z + _const_fm.z);
      double fzz = wn < 0. ? val(f,i,0,1) - val(f,i,0,0) : val(f,i,0,0) - val(f,i,0,-1);
      f2 -= dt*wn*fzz/(2.*Delta);
    }


    val(flux.x,0,0,0) = f2*val(uf.x,0,0,0);
  }}end_is_face_x()
#line 28
is_face_y(){ {







    double un = dt*val(uf.y,0,0,0)/(_const_fm.y*Delta + 0.), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = val(f,0,i,0) + (_const_src + _const_src)*dt/4. + s*(1. - s*un)*val(g.y,0,i,0)*Delta/2.;





    if (_const_fm.z && _const_fm.z) {
      double vn = (val(uf.z,0,i,0) + val(uf.z,0,i,1))/(_const_fm.z + _const_fm.z);
      double fyy = vn < 0. ? val(f,0,i,1) - val(f,0,i,0) : val(f,0,i,0) - val(f,0,i,-1);
      f2 -= dt*vn*fyy/(2.*Delta);
    }


    if (_const_fm.x && _const_fm.x) {
      double wn = (val(uf.x,0,i,0) + val(uf.x,1,i,0))/(_const_fm.x + _const_fm.x);
      double fzz = wn < 0. ? val(f,1,i,0) - val(f,0,i,0) : val(f,0,i,0) - val(f,-1,i,0);
      f2 -= dt*wn*fzz/(2.*Delta);
    }


    val(flux.y,0,0,0) = f2*val(uf.y,0,0,0);
  }}end_is_face_y()
#line 28
is_face_z(){ {







    double un = dt*val(uf.z,0,0,0)/(_const_fm.z*Delta + 0.), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = val(f,0,0,i) + (_const_src + _const_src)*dt/4. + s*(1. - s*un)*val(g.z,0,0,i)*Delta/2.;





    if (_const_fm.x && _const_fm.x) {
      double vn = (val(uf.x,0,0,i) + val(uf.x,1,0,i))/(_const_fm.x + _const_fm.x);
      double fyy = vn < 0. ? val(f,1,0,i) - val(f,0,0,i) : val(f,0,0,i) - val(f,-1,0,i);
      f2 -= dt*vn*fyy/(2.*Delta);
    }


    if (_const_fm.y && _const_fm.y) {
      double wn = (val(uf.y,0,0,i) + val(uf.y,0,1,i))/(_const_fm.y + _const_fm.y);
      double fzz = wn < 0. ? val(f,0,1,i) - val(f,0,0,i) : val(f,0,0,i) - val(f,0,-1,i);
      f2 -= dt*wn*fzz/(2.*Delta);
    }


    val(flux.z,0,0,0) = f2*val(uf.z,0,0,0);
  }}end_is_face_z()}end_foreach_face_generic();}}delete((scalar*)((vector[]){g,{{-1},{-1},{-1}}}));
}






void advection (scalar * tracers, vector u, double dt,
  scalar * src)
{




  scalar * psrc = src;
  if (!src)
    {scalar*_i=(scalar*)( tracers);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {
      const scalar zero = new_const_scalar("zero",12, 0.);
      src = list_append (src, zero);
    }}}
  if (!(list_len (tracers) == list_len (src))) qassert ("/home/lisergey/basilisk/src/bcg.h", 80, "list_len (tracers) == list_len (src)");

  scalar f, source;
  {scalar*_i0=src;scalar*_i1= tracers;if(_i0)for(source=*_i0,f=*_i1;_i0->i>= 0;source=*++_i0,f=*++_i1){ {
    vector  flux=new_face_vector("flux");
    tracer_fluxes (f, u, flux, dt, source);

    if(!is_constant(cm)){

    
#line 87
foreach_stencil(1,{(NonLocal[]){{"cm","scalar",(void *)&cm,NULL,0},{"flux","vector",(void *)&flux,NULL,0},{"dt","double",(void *)&dt,NULL,0},{"f","scalar",(void *)&f,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n      \n// #line 88 \"/home/lisergey/basilisk/src/bcg.h\"\n{\n        val_out_(f,0,0,0) += dt*(val(flux.x,0,0,0) - val(flux.x,1,0,0))/(Delta*val(cm,0,0,0));\n        \n// #line 89\nval_out_(f,0,0,0) += dt*(val(flux.y,0,0,0) - val(flux.y,0,1,0))/(Delta*val(cm,0,0,0));\n        \n// #line 89\nval_out_(f,0,0,0) += dt*(val(flux.z,0,0,0) - val(flux.z,0,0,1))/(Delta*val(cm,0,0,0));}"})
      {
        {_stencil_val(flux.x,0,0,0); _stencil_val(flux.x,1,0,0);_stencil_val(cm,0,0,0);_stencil_val_r(f,0,0,0);   }
        
#line 89
{_stencil_val(flux.y,0,0,0); _stencil_val(flux.y,0,1,0);_stencil_val(cm,0,0,0);_stencil_val_r(f,0,0,0);   }
        
#line 89
{_stencil_val(flux.z,0,0,0); _stencil_val(flux.z,0,0,1);_stencil_val(cm,0,0,0);_stencil_val_r(f,0,0,0);   }}end_foreach_stencil();{
#line 87
foreach()
      {
        val(f,0,0,0) += dt*(val(flux.x,0,0,0) - val(flux.x,1,0,0))/(Delta*val(cm,0,0,0));
        
#line 89
val(f,0,0,0) += dt*(val(flux.y,0,0,0) - val(flux.y,0,1,0))/(Delta*val(cm,0,0,0));
        
#line 89
val(f,0,0,0) += dt*(val(flux.z,0,0,0) - val(flux.z,0,0,1))/(Delta*val(cm,0,0,0));}end_foreach();}}else {double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);

    
#line 87
foreach_stencil(1,{(NonLocal[]){{"_const_cm","double",(void *)&_const_cm,NULL,0},{"flux","vector",(void *)&flux,NULL,0},{"dt","double",(void *)&dt,NULL,0},{"f","scalar",(void *)&f,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n      \n// #line 88 \"/home/lisergey/basilisk/src/bcg.h\"\n{\n        val_out_(f,0,0,0) += dt*(val(flux.x,0,0,0) - val(flux.x,1,0,0))/(Delta*_const_cm);\n        \n// #line 89\nval_out_(f,0,0,0) += dt*(val(flux.y,0,0,0) - val(flux.y,0,1,0))/(Delta*_const_cm);\n        \n// #line 89\nval_out_(f,0,0,0) += dt*(val(flux.z,0,0,0) - val(flux.z,0,0,1))/(Delta*_const_cm);}"})
      {
        {_stencil_val(flux.x,0,0,0); _stencil_val(flux.x,1,0,0);;_stencil_val_r(f,0,0,0);   }
        
#line 89
{_stencil_val(flux.y,0,0,0); _stencil_val(flux.y,0,1,0);;_stencil_val_r(f,0,0,0);   }
        
#line 89
{_stencil_val(flux.z,0,0,0); _stencil_val(flux.z,0,0,1);;_stencil_val_r(f,0,0,0);   }}end_foreach_stencil();

    {
#line 87
foreach()
      {
        val(f,0,0,0) += dt*(val(flux.x,0,0,0) - val(flux.x,1,0,0))/(Delta*_const_cm);
        
#line 89
val(f,0,0,0) += dt*(val(flux.y,0,0,0) - val(flux.y,0,1,0))/(Delta*_const_cm);
        
#line 89
val(f,0,0,0) += dt*(val(flux.z,0,0,0) - val(flux.z,0,0,1))/(Delta*_const_cm);}end_foreach();}}delete((scalar*)((vector[]){flux,{{-1},{-1},{-1}}}));



  }}}

  if (!psrc)
    pfree (src,__func__,__FILE__,__LINE__);
}
#line 30 "/home/lisergey/basilisk/src/navier-stokes/centered.h"



#line 1 "./viscosity.h"
#line 1 "/home/lisergey/basilisk/src/viscosity.h"
#line 50 "/home/lisergey/basilisk/src/viscosity.h"
#line 1 "./poisson.h"
#line 1 "/home/lisergey/basilisk/src/poisson.h"
#line 32 "/home/lisergey/basilisk/src/poisson.h"
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
  
     val(s,0,0,0) = bilinear (point, s);}}end_foreach_level();}





    boundary_level (da, l);
    for (int i = 0; i < nrelax; i++) {
      relax (da, res, l, data);
      boundary_level (da, l);
    }
  }




  foreach_stencil(1,{(NonLocal[]){{"da","scalar",(void *)da,NULL,1},{"a","scalar",(void *)a,NULL,1},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},""," \n// #line 84 \"/home/lisergey/basilisk/src/poisson.h\"\n{\n    scalar s, ds;\n    {forin2 (s, ds , a, da)\n     \n val_out_(s,0,0,0) += val(ds,0,0,0); endforin2()}\n  }"}) {
    scalar s, ds;
    {scalar*_i0= da;scalar*_i1= a;if(_i0)for(ds=*_i0,s=*_i1;_i0->i>= 0;ds=*++_i0,s=*++_i1){
     
 { _stencil_val(ds,0,0,0);_stencil_val_r(s,0,0,0); }}}
  }end_foreach_stencil();




  {
#line 84
foreach() {
    scalar s, ds;
    {scalar*_i0= da;scalar*_i1= a;if(_i0)for(ds=*_i0,s=*_i1;_i0->i>= 0;ds=*++_i0,s=*++_i1){
     
 val(s,0,0,0) += val(ds,0,0,0);}}
  }end_foreach();}
}
#line 102 "/home/lisergey/basilisk/src/poisson.h"
int NITERMAX = 100, NITERMIN = 1;
double TOLERANCE = 1e-3;




typedef struct {
  int i;
  double resb, resa;
  double sum;
  int nrelax;
  int minlevel;
} mgstats;
#line 125 "/home/lisergey/basilisk/src/poisson.h"
mgstats mg_solve (scalar * a, scalar * b,
    double (* residual) (scalar * a, scalar * b, scalar * res,
           void * data),
    void (* relax) (scalar * da, scalar * res, int depth,
      void * data),
    void * data,
    int nrelax,
    scalar * res,
    int minlevel,
    double tolerance)
{





  scalar * da = list_clone (a), * pres = res;
  if (!res)
    res = list_clone (b);






  for (int b = 0; b < nboundary; b++)
    {scalar*_i=(scalar*)( da);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      _attribute[s.i].boundary[b] = _attribute[s.i].boundary_homogeneous[b];}}




  mgstats s = {0};
  double sum = 0.;
  scalar rhs = b[0];
  foreach_stencil (1,{(NonLocal[]){{"rhs","scalar",(void *)&rhs,NULL,0},{"sum","double",(void *)&sum,NULL,0,'+'},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n    \n// #line 161 \"/home/lisergey/basilisk/src/poisson.h\"\nsum += val(rhs,0,0,0);"})
    { _stencil_val(rhs,0,0,0); }end_foreach_stencil();
  
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel reduction(+:sum)){
#line 160
foreach ()
    sum += val(rhs,0,0,0);end_foreach();mpi_all_reduce_array(&sum,double,MPI_SUM,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
  
#line 162
s.sum = sum;
  s.nrelax = nrelax > 0 ? nrelax : 4;




  double resb;
  resb = s.resb = s.resa = (* residual) (a, b, res, data);






  for (s.i = 0;
       s.i < NITERMAX && (s.i < NITERMIN || s.resa > tolerance);
       s.i++) {
    mg_cycle (a, res, da, relax, data,
       s.nrelax,
       minlevel,
       grid->maxdepth);
    s.resa = (* residual) (a, b, res, data);
#line 192 "/home/lisergey/basilisk/src/poisson.h"
    if (s.resa > tolerance) {
      if (resb/s.resa < 1.2 && s.nrelax < 100)
 s.nrelax++;
      else if (resb/s.resa > 10 && s.nrelax > 2)
 s.nrelax--;
    }







    resb = s.resa;
  }
  s.minlevel = minlevel;




  if (s.resa > tolerance) {
    scalar v = a[0];
    fprintf (ferr,
      "WARNING: convergence for %s not reached after %d iterations\n"
      "  res: %g sum: %g nrelax: %d tolerance: %g\n", _attribute[v.i].name,
      s.i, s.resa, s.sum, s.nrelax, tolerance), fflush (ferr);
  }




  if (!pres)
    delete (res), pfree (res,__func__,__FILE__,__LINE__);
  delete (da), pfree (da,__func__,__FILE__,__LINE__);

  return s;
}
#line 251 "/home/lisergey/basilisk/src/poisson.h"
struct Poisson {
  scalar a, b;
          vector alpha;
          scalar lambda;
  double tolerance;
  int nrelax, minlevel;
  scalar * res;



};





static void relax (scalar * al, scalar * bl, int l, void * data)
{
  scalar a = al[0], b = bl[0];
  struct Poisson * p = (struct Poisson *) data;
          vector alpha = p->alpha;
          scalar lambda = p->lambda;
#line 289 "/home/lisergey/basilisk/src/poisson.h"
  scalar c = a;






  if(!is_constant(lambda) && !is_constant(alpha.x)){{foreach_level_or_leaf (l) {
    double n = - sq(Delta)*val(b,0,0,0), d = - val(lambda,0,0,0)*sq(Delta);
     {
      n += val(alpha.x,1,0,0)*val(a,1,0,0) + val(alpha.x,0,0,0)*val(a,-1,0,0);
      d += val(alpha.x,1,0,0) + val(alpha.x,0,0,0);
    } 
#line 298
{
      n += val(alpha.y,0,1,0)*val(a,0,1,0) + val(alpha.y,0,0,0)*val(a,0,-1,0);
      d += val(alpha.y,0,1,0) + val(alpha.y,0,0,0);
    } 
#line 298
{
      n += val(alpha.z,0,0,1)*val(a,0,0,1) + val(alpha.z,0,0,0)*val(a,0,0,-1);
      d += val(alpha.z,0,0,1) + val(alpha.z,0,0,0);
    }
#line 312 "/home/lisergey/basilisk/src/poisson.h"
      val(c,0,0,0) = n/d;
  }end_foreach_level_or_leaf();}}else if(is_constant(lambda) && !is_constant(alpha.x)){double _const_lambda=_constant[lambda.i-_NVARMAX];NOT_UNUSED(_const_lambda);






  {
#line 296
foreach_level_or_leaf (l) {
    double n = - sq(Delta)*val(b,0,0,0), d = - _const_lambda*sq(Delta);
     {
      n += val(alpha.x,1,0,0)*val(a,1,0,0) + val(alpha.x,0,0,0)*val(a,-1,0,0);
      d += val(alpha.x,1,0,0) + val(alpha.x,0,0,0);
    } 
#line 298
{
      n += val(alpha.y,0,1,0)*val(a,0,1,0) + val(alpha.y,0,0,0)*val(a,0,-1,0);
      d += val(alpha.y,0,1,0) + val(alpha.y,0,0,0);
    } 
#line 298
{
      n += val(alpha.z,0,0,1)*val(a,0,0,1) + val(alpha.z,0,0,0)*val(a,0,0,-1);
      d += val(alpha.z,0,0,1) + val(alpha.z,0,0,0);
    }
#line 312 "/home/lisergey/basilisk/src/poisson.h"
      val(c,0,0,0) = n/d;
  }end_foreach_level_or_leaf();}}else if(!is_constant(lambda) && is_constant(alpha.x)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX],_constant[alpha.z.i-_NVARMAX]};NOT_UNUSED(_const_alpha);






  {
#line 296
foreach_level_or_leaf (l) {
    double n = - sq(Delta)*val(b,0,0,0), d = - val(lambda,0,0,0)*sq(Delta);
     {
      n += _const_alpha.x*val(a,1,0,0) + _const_alpha.x*val(a,-1,0,0);
      d += _const_alpha.x + _const_alpha.x;
    } 
#line 298
{
      n += _const_alpha.y*val(a,0,1,0) + _const_alpha.y*val(a,0,-1,0);
      d += _const_alpha.y + _const_alpha.y;
    } 
#line 298
{
      n += _const_alpha.z*val(a,0,0,1) + _const_alpha.z*val(a,0,0,-1);
      d += _const_alpha.z + _const_alpha.z;
    }
#line 312 "/home/lisergey/basilisk/src/poisson.h"
      val(c,0,0,0) = n/d;
  }end_foreach_level_or_leaf();}}else {double _const_lambda=_constant[lambda.i-_NVARMAX];NOT_UNUSED(_const_lambda);_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX],_constant[alpha.z.i-_NVARMAX]};NOT_UNUSED(_const_alpha);






  {
#line 296
foreach_level_or_leaf (l) {
    double n = - sq(Delta)*val(b,0,0,0), d = - _const_lambda*sq(Delta);
     {
      n += _const_alpha.x*val(a,1,0,0) + _const_alpha.x*val(a,-1,0,0);
      d += _const_alpha.x + _const_alpha.x;
    } 
#line 298
{
      n += _const_alpha.y*val(a,0,1,0) + _const_alpha.y*val(a,0,-1,0);
      d += _const_alpha.y + _const_alpha.y;
    } 
#line 298
{
      n += _const_alpha.z*val(a,0,0,1) + _const_alpha.z*val(a,0,0,-1);
      d += _const_alpha.z + _const_alpha.z;
    }
#line 312 "/home/lisergey/basilisk/src/poisson.h"
      val(c,0,0,0) = n/d;
  }end_foreach_level_or_leaf();}}
#line 331 "/home/lisergey/basilisk/src/poisson.h"
}






static double residual (scalar * al, scalar * bl, scalar * resl, void * data)
{
  scalar a = al[0], b = bl[0], res = resl[0];
  struct Poisson * p = (struct Poisson *) data;
          vector alpha = p->alpha;
          scalar lambda = p->lambda;
  double maxres = 0.;


  vector  g=new_face_vector("g");
  if(!is_constant(alpha.x)){
  
#line 348
foreach_face_stencil(1,{(NonLocal[]){{"a","scalar",(void *)&a,NULL,0},{"alpha","vector",(void *)&alpha,NULL,0},{"g","vector",(void *)&g,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n// #line 348 \"/home/lisergey/basilisk/src/poisson.h\"\n{is_face_x(){\n    val_out_(g.x,0,0,0) = val(alpha.x,0,0,0)*((val(a,0,0,0) - val(a,0 -1,0,0))/Delta);}end_is_face_x()\n// #line 348\nis_face_y(){\n    val_out_(g.y,0,0,0) = val(alpha.y,0,0,0)*((val(a,0,0,0) - val(a,0,0 -1,0))/Delta);}end_is_face_y()\n// #line 348\nis_face_z(){\n    val_out_(g.z,0,0,0) = val(alpha.z,0,0,0)*((val(a,0,0,0) - val(a,0,0,0 -1))/Delta);}end_is_face_z()}"}){_stencil_is_face_x(){
    { _stencil_val(alpha.x,0,0,0);_stencil_val(a,0,0,0); _stencil_val(a,0 -1,0,0);_stencil_val_a(g.x,0,0,0);  }}end__stencil_is_face_x()
#line 348
_stencil_is_face_y(){
    { _stencil_val(alpha.y,0,0,0);_stencil_val(a,0,0,0); _stencil_val(a,0,0 -1,0);_stencil_val_a(g.y,0,0,0);  }}end__stencil_is_face_y()
#line 348
_stencil_is_face_z(){
    { _stencil_val(alpha.z,0,0,0);_stencil_val(a,0,0,0); _stencil_val(a,0,0,0 -1);_stencil_val_a(g.z,0,0,0);  }}end__stencil_is_face_z()}end_foreach_face_stencil();{
#line 348
foreach_face_generic(){is_face_x(){
    val(g.x,0,0,0) = val(alpha.x,0,0,0)*((val(a,0,0,0) - val(a,0 -1,0,0))/Delta);}end_is_face_x()
#line 348
is_face_y(){
    val(g.y,0,0,0) = val(alpha.y,0,0,0)*((val(a,0,0,0) - val(a,0,0 -1,0))/Delta);}end_is_face_y()
#line 348
is_face_z(){
    val(g.z,0,0,0) = val(alpha.z,0,0,0)*((val(a,0,0,0) - val(a,0,0,0 -1))/Delta);}end_is_face_z()}end_foreach_face_generic();}}else {_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX],_constant[alpha.z.i-_NVARMAX]};NOT_UNUSED(_const_alpha);
  
#line 348
foreach_face_stencil(1,{(NonLocal[]){{"a","scalar",(void *)&a,NULL,0},{"_const_alpha","_coord",(void *)&_const_alpha,NULL,0},{"g","vector",(void *)&g,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n// #line 348 \"/home/lisergey/basilisk/src/poisson.h\"\n{is_face_x(){\n    val_out_(g.x,0,0,0) = _const_alpha.x*((val(a,0,0,0) - val(a,0 -1,0,0))/Delta);}end_is_face_x()\n// #line 348\nis_face_y(){\n    val_out_(g.y,0,0,0) = _const_alpha.y*((val(a,0,0,0) - val(a,0,0 -1,0))/Delta);}end_is_face_y()\n// #line 348\nis_face_z(){\n    val_out_(g.z,0,0,0) = _const_alpha.z*((val(a,0,0,0) - val(a,0,0,0 -1))/Delta);}end_is_face_z()}"}){_stencil_is_face_x(){
    {;_stencil_val(a,0,0,0); _stencil_val(a,0 -1,0,0);_stencil_val_a(g.x,0,0,0);  }}end__stencil_is_face_x()
#line 348
_stencil_is_face_y(){
    {;_stencil_val(a,0,0,0); _stencil_val(a,0,0 -1,0);_stencil_val_a(g.y,0,0,0);  }}end__stencil_is_face_y()
#line 348
_stencil_is_face_z(){
    {;_stencil_val(a,0,0,0); _stencil_val(a,0,0,0 -1);_stencil_val_a(g.z,0,0,0);  }}end__stencil_is_face_z()}end_foreach_face_stencil();
  {
#line 348
foreach_face_generic(){is_face_x(){
    val(g.x,0,0,0) = _const_alpha.x*((val(a,0,0,0) - val(a,0 -1,0,0))/Delta);}end_is_face_x()
#line 348
is_face_y(){
    val(g.y,0,0,0) = _const_alpha.y*((val(a,0,0,0) - val(a,0,0 -1,0))/Delta);}end_is_face_y()
#line 348
is_face_z(){
    val(g.z,0,0,0) = _const_alpha.z*((val(a,0,0,0) - val(a,0,0,0 -1))/Delta);}end_is_face_z()}end_foreach_face_generic();}}
  if(!is_constant(lambda)){
  
#line 350
foreach_stencil (1,{(NonLocal[]){{"maxres","double",(void *)&maxres,NULL,0,'M'},{"g","vector",(void *)&g,NULL,0},{"a","scalar",(void *)&a,NULL,0},{"lambda","scalar",(void *)&lambda,NULL,0},{"b","scalar",(void *)&b,NULL,0},{"res","scalar",(void *)&res,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},""," \n// #line 350 \"/home/lisergey/basilisk/src/poisson.h\"\n{\n    val_out_(res,0,0,0) = val(b,0,0,0) - val(lambda,0,0,0)*val(a,0,0,0);\n    \n      val_out_(res,0,0,0) -= (val(g.x,1,0,0) - val(g.x,0,0,0))/Delta;\n      \n// #line 353\nval_out_(res,0,0,0) -= (val(g.y,0,1,0) - val(g.y,0,0,0))/Delta;\n      \n// #line 353\nval_out_(res,0,0,0) -= (val(g.z,0,0,1) - val(g.z,0,0,0))/Delta;\n\n\n\n\n\n\n    if (fabs (val(res,0,0,0)) > maxres)\n      maxres = fabs (val(res,0,0,0));\n  }"}) { 
_stencil_val(b,0,0,0); _stencil_val(lambda,0,0,0);_stencil_val(a,0,0,0);
    
#line 351
_stencil_val_a(res,0,0,0);  
    
      {_stencil_val(g.x,1,0,0); _stencil_val(g.x,0,0,0);_stencil_val_r(res,0,0,0);   }
      
#line 353
{_stencil_val(g.y,0,1,0); _stencil_val(g.y,0,0,0);_stencil_val_r(res,0,0,0);   }
      
#line 353
{_stencil_val(g.z,0,0,1); _stencil_val(g.z,0,0,0);_stencil_val_r(res,0,0,0);   }






_stencil_val(res,0,0,0);
      {_stencil_val(res,0,0,0);   }






        
  
#line 362
}end_foreach_stencil();
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel reduction(max:maxres)){
#line 350
foreach () {
    val(res,0,0,0) = val(b,0,0,0) - val(lambda,0,0,0)*val(a,0,0,0);
    
      val(res,0,0,0) -= (val(g.x,1,0,0) - val(g.x,0,0,0))/Delta;
      
#line 353
val(res,0,0,0) -= (val(g.y,0,1,0) - val(g.y,0,0,0))/Delta;
      
#line 353
val(res,0,0,0) -= (val(g.z,0,0,1) - val(g.z,0,0,0))/Delta;






    if (fabs (val(res,0,0,0)) > maxres)
      maxres = fabs (val(res,0,0,0));
  }end_foreach();mpi_all_reduce_array(&maxres,double,MPI_MAX,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 362
}else {double _const_lambda=_constant[lambda.i-_NVARMAX];NOT_UNUSED(_const_lambda);
  
#line 350
foreach_stencil (1,{(NonLocal[]){{"maxres","double",(void *)&maxres,NULL,0,'M'},{"g","vector",(void *)&g,NULL,0},{"a","scalar",(void *)&a,NULL,0},{"_const_lambda","double",(void *)&_const_lambda,NULL,0},{"b","scalar",(void *)&b,NULL,0},{"res","scalar",(void *)&res,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},""," \n// #line 350 \"/home/lisergey/basilisk/src/poisson.h\"\n{\n    val_out_(res,0,0,0) = val(b,0,0,0) - _const_lambda*val(a,0,0,0);\n    \n      val_out_(res,0,0,0) -= (val(g.x,1,0,0) - val(g.x,0,0,0))/Delta;\n      \n// #line 353\nval_out_(res,0,0,0) -= (val(g.y,0,1,0) - val(g.y,0,0,0))/Delta;\n      \n// #line 353\nval_out_(res,0,0,0) -= (val(g.z,0,0,1) - val(g.z,0,0,0))/Delta;\n\n\n\n\n\n\n    if (fabs (val(res,0,0,0)) > maxres)\n      maxres = fabs (val(res,0,0,0));\n  }"}) { 
_stencil_val(b,0,0,0);;_stencil_val(a,0,0,0);
    
#line 351
_stencil_val_a(res,0,0,0);  
    
      {_stencil_val(g.x,1,0,0); _stencil_val(g.x,0,0,0);_stencil_val_r(res,0,0,0);   }
      
#line 353
{_stencil_val(g.y,0,1,0); _stencil_val(g.y,0,0,0);_stencil_val_r(res,0,0,0);   }
      
#line 353
{_stencil_val(g.z,0,0,1); _stencil_val(g.z,0,0,0);_stencil_val_r(res,0,0,0);   }






_stencil_val(res,0,0,0);
      {_stencil_val(res,0,0,0);   }






        
  
#line 362
}end_foreach_stencil();
  
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel reduction(max:maxres)){
#line 350
foreach () {
    val(res,0,0,0) = val(b,0,0,0) - _const_lambda*val(a,0,0,0);
    
      val(res,0,0,0) -= (val(g.x,1,0,0) - val(g.x,0,0,0))/Delta;
      
#line 353
val(res,0,0,0) -= (val(g.y,0,1,0) - val(g.y,0,0,0))/Delta;
      
#line 353
val(res,0,0,0) -= (val(g.z,0,0,1) - val(g.z,0,0,0))/Delta;






    if (fabs (val(res,0,0,0)) > maxres)
      maxres = fabs (val(res,0,0,0));
  }end_foreach();mpi_all_reduce_array(&maxres,double,MPI_MAX,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 362
}
#line 380 "/home/lisergey/basilisk/src/poisson.h"
  {delete((scalar*)((vector[]){g,{{-1},{-1},{-1}}}));return maxres;}delete((scalar*)((vector[]){g,{{-1},{-1},{-1}}}));
}
#line 392 "/home/lisergey/basilisk/src/poisson.h"
mgstats poisson (scalar a, scalar b,
           vector alpha,
           scalar lambda,
   double tolerance,
   int nrelax,
   int minlevel,
   scalar * res,
   double (* flux) (Point, scalar, vector, double *))
{






  if (alpha.x.i < 0)
    alpha = unityf;
  if (lambda.i < 0) {
    const scalar zeroc = new_const_scalar("zeroc",13, 0.);
    lambda = zeroc;
  }




  restriction (((scalar[]){alpha.x,alpha.y,alpha.z,lambda,{-1}}));





  double defaultol = TOLERANCE;
  if (tolerance)
    TOLERANCE = tolerance;

  struct Poisson p = {a, b, alpha, lambda, tolerance, nrelax, minlevel, res };






  mgstats s = mg_solve ((
#line 125
scalar *
#line 434
)((scalar[]){a,{-1}}),( 
#line 125
scalar *
#line 434
)((scalar[]) {b,{-1}}), residual, relax, &p
,
   
#line 435
nrelax, res, max(1, minlevel)
#line 133
, 
TOLERANCE
#line 435
);




  if (tolerance)
    TOLERANCE = defaultol;

  return s;
}
#line 463 "/home/lisergey/basilisk/src/poisson.h"
     
mgstats project (vector uf, scalar p,
           vector alpha,
   double dt,
   int nrelax)
{tracing("project","/home/lisergey/basilisk/src/poisson.h",464);






  scalar  div=new_scalar("div");
  foreach_stencil(1,{(NonLocal[]){{"dt","double",(void *)&dt,NULL,0},{"uf","vector",(void *)&uf,NULL,0},{"div","scalar",(void *)&div,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},""," \n// #line 476 \"/home/lisergey/basilisk/src/poisson.h\"\n{\n    val_out_(div,0,0,0) = 0.;\n    \n      val_out_(div,0,0,0) += val(uf.x,1,0,0) - val(uf.x,0,0,0);\n      \n// #line 479\nval_out_(div,0,0,0) += val(uf.y,0,1,0) - val(uf.y,0,0,0);\n      \n// #line 479\nval_out_(div,0,0,0) += val(uf.z,0,0,1) - val(uf.z,0,0,0);\n    val_out_(div,0,0,0) /= dt*Delta;\n  }"}) {
    _stencil_val_a(div,0,0,0);  
    
      { _stencil_val(uf.x,1,0,0); _stencil_val(uf.x,0,0,0);_stencil_val_r(div,0,0,0);  }
      
#line 479
{ _stencil_val(uf.y,0,1,0); _stencil_val(uf.y,0,0,0);_stencil_val_r(div,0,0,0);  }
      
#line 479
{ _stencil_val(uf.z,0,0,1); _stencil_val(uf.z,0,0,0);_stencil_val_r(div,0,0,0);  }
    _stencil_val_r(div,0,0,0);  
  }end_foreach_stencil();
  {
#line 476
foreach() {
    val(div,0,0,0) = 0.;
    
      val(div,0,0,0) += val(uf.x,1,0,0) - val(uf.x,0,0,0);
      
#line 479
val(div,0,0,0) += val(uf.y,0,1,0) - val(uf.y,0,0,0);
      
#line 479
val(div,0,0,0) += val(uf.z,0,0,1) - val(uf.z,0,0,0);
    val(div,0,0,0) /= dt*Delta;
  }end_foreach();}
#line 492 "/home/lisergey/basilisk/src/poisson.h"
  mgstats mgp = poisson (p, div, alpha
#line 393
,
( scalar) {-1}
#line 493
, TOLERANCE/sq(dt), nrelax
#line 396
, 
0, 
NULL, 
NULL
#line 493
);




  if(!is_constant(alpha.x)){




  
#line 498
foreach_face_stencil(1,{(NonLocal[]){{"p","scalar",(void *)&p,NULL,0},{"alpha","vector",(void *)&alpha,NULL,0},{"dt","double",(void *)&dt,NULL,0},{"uf","vector",(void *)&uf,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n// #line 498 \"/home/lisergey/basilisk/src/poisson.h\"\n{is_face_x(){\n    val_out_(uf.x,0,0,0) -= dt*val(alpha.x,0,0,0)*((val(p,0,0,0) - val(p,0 -1,0,0))/Delta);}end_is_face_x()\n// #line 498\nis_face_y(){\n    val_out_(uf.y,0,0,0) -= dt*val(alpha.y,0,0,0)*((val(p,0,0,0) - val(p,0,0 -1,0))/Delta);}end_is_face_y()\n// #line 498\nis_face_z(){\n    val_out_(uf.z,0,0,0) -= dt*val(alpha.z,0,0,0)*((val(p,0,0,0) - val(p,0,0,0 -1))/Delta);}end_is_face_z()}"}){_stencil_is_face_x(){
    {_stencil_val(alpha.x,0,0,0);_stencil_val(p,0,0,0); _stencil_val(p,0 -1,0,0);_stencil_val_r(uf.x,0,0,0);   }}end__stencil_is_face_x()
#line 498
_stencil_is_face_y(){
    {_stencil_val(alpha.y,0,0,0);_stencil_val(p,0,0,0); _stencil_val(p,0,0 -1,0);_stencil_val_r(uf.y,0,0,0);   }}end__stencil_is_face_y()
#line 498
_stencil_is_face_z(){
    {_stencil_val(alpha.z,0,0,0);_stencil_val(p,0,0,0); _stencil_val(p,0,0,0 -1);_stencil_val_r(uf.z,0,0,0);   }}end__stencil_is_face_z()}end_foreach_face_stencil();{
#line 498
foreach_face_generic(){is_face_x(){
    val(uf.x,0,0,0) -= dt*val(alpha.x,0,0,0)*((val(p,0,0,0) - val(p,0 -1,0,0))/Delta);}end_is_face_x()
#line 498
is_face_y(){
    val(uf.y,0,0,0) -= dt*val(alpha.y,0,0,0)*((val(p,0,0,0) - val(p,0,0 -1,0))/Delta);}end_is_face_y()
#line 498
is_face_z(){
    val(uf.z,0,0,0) -= dt*val(alpha.z,0,0,0)*((val(p,0,0,0) - val(p,0,0,0 -1))/Delta);}end_is_face_z()}end_foreach_face_generic();}}else {_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX],_constant[alpha.z.i-_NVARMAX]};NOT_UNUSED(_const_alpha);




  
#line 498
foreach_face_stencil(1,{(NonLocal[]){{"p","scalar",(void *)&p,NULL,0},{"_const_alpha","_coord",(void *)&_const_alpha,NULL,0},{"dt","double",(void *)&dt,NULL,0},{"uf","vector",(void *)&uf,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n// #line 498 \"/home/lisergey/basilisk/src/poisson.h\"\n{is_face_x(){\n    val_out_(uf.x,0,0,0) -= dt*_const_alpha.x*((val(p,0,0,0) - val(p,0 -1,0,0))/Delta);}end_is_face_x()\n// #line 498\nis_face_y(){\n    val_out_(uf.y,0,0,0) -= dt*_const_alpha.y*((val(p,0,0,0) - val(p,0,0 -1,0))/Delta);}end_is_face_y()\n// #line 498\nis_face_z(){\n    val_out_(uf.z,0,0,0) -= dt*_const_alpha.z*((val(p,0,0,0) - val(p,0,0,0 -1))/Delta);}end_is_face_z()}"}){_stencil_is_face_x(){
    {;_stencil_val(p,0,0,0); _stencil_val(p,0 -1,0,0);_stencil_val_r(uf.x,0,0,0);   }}end__stencil_is_face_x()
#line 498
_stencil_is_face_y(){
    {;_stencil_val(p,0,0,0); _stencil_val(p,0,0 -1,0);_stencil_val_r(uf.y,0,0,0);   }}end__stencil_is_face_y()
#line 498
_stencil_is_face_z(){
    {;_stencil_val(p,0,0,0); _stencil_val(p,0,0,0 -1);_stencil_val_r(uf.z,0,0,0);   }}end__stencil_is_face_z()}end_foreach_face_stencil();




  {
#line 498
foreach_face_generic(){is_face_x(){
    val(uf.x,0,0,0) -= dt*_const_alpha.x*((val(p,0,0,0) - val(p,0 -1,0,0))/Delta);}end_is_face_x()
#line 498
is_face_y(){
    val(uf.y,0,0,0) -= dt*_const_alpha.y*((val(p,0,0,0) - val(p,0,0 -1,0))/Delta);}end_is_face_y()
#line 498
is_face_z(){
    val(uf.z,0,0,0) -= dt*_const_alpha.z*((val(p,0,0,0) - val(p,0,0,0 -1))/Delta);}end_is_face_z()}end_foreach_face_generic();}}

  {delete((scalar*)((scalar[]){div,{-1}}));{end_tracing("project","/home/lisergey/basilisk/src/poisson.h",501);return mgp;}}delete((scalar*)((scalar[]){div,{-1}}));
end_tracing("project","/home/lisergey/basilisk/src/poisson.h",502);}
#line 51 "/home/lisergey/basilisk/src/viscosity.h"

struct Viscosity {
  vector mu;
  scalar rho;
  double dt;
};
#line 135 "/home/lisergey/basilisk/src/viscosity.h"
static void relax_viscosity (scalar * a, scalar * b, int l, void * data)
{
  struct Viscosity * p = (struct Viscosity *) data;
          vector mu = p->mu;
          scalar rho = p->rho;
  double dt = p->dt;
  vector u = (*((vector *)&(a[0]))), r = (*((vector *)&(b[0])));




  vector w = u;


  if(!is_constant(rho) && !is_constant(mu.x)){{foreach_level_or_leaf (l) {
    
      val(w.x,0,0,0) = (dt/val(rho,0,0,0)*(2.*val(mu.x,1,0,0)*val(u.x,1,0,0) + 2.*val(mu.x,0,0,0)*val(u.x,-1,0,0)

      + val(mu.y,0,1,0)*(val(u.x,0,1,0) +
     (val(u.y,1,0,0) + val(u.y,1,1,0))/4. -
     (val(u.y,-1,0,0) + val(u.y,-1,1,0))/4.)
      - val(mu.y,0,0,0)*(- val(u.x,0,-1,0) +
         (val(u.y,1,-1,0) + val(u.y,1,0,0))/4. -
         (val(u.y,-1,-1,0) + val(u.y,-1,0,0))/4.)


      + val(mu.z,0,0,1)*(val(u.x,0,0,1) +
       (val(u.z,1,0,0) + val(u.z,1,0,1))/4. -
       (val(u.z,-1,0,0) + val(u.z,-1,0,1))/4.)
      - val(mu.z,0,0,0)*(- val(u.x,0,0,-1) +
         (val(u.z,1,0,-1) + val(u.z,1,0,0))/4. -
         (val(u.z,-1,0,-1) + val(u.z,-1,0,0))/4.)

      ) + val(r.x,0,0,0)*sq(Delta))/
    (sq(Delta)*((coord){1.,1.,1.}).x + dt/val(rho,0,0,0)*(2.*val(mu.x,1,0,0) + 2.*val(mu.x,0,0,0)

          + val(mu.y,0,1,0) + val(mu.y,0,0,0)


          + val(mu.z,0,0,1) + val(mu.z,0,0,0)

        ));
      
#line 151
val(w.y,0,0,0) = (dt/val(rho,0,0,0)*(2.*val(mu.y,0,1,0)*val(u.y,0,1,0) + 2.*val(mu.y,0,0,0)*val(u.y,0,-1,0)

      + val(mu.z,0,0,1)*(val(u.y,0,0,1) +
     (val(u.z,0,1,0) + val(u.z,0,1,1))/4. -
     (val(u.z,0,-1,0) + val(u.z,0,-1,1))/4.)
      - val(mu.z,0,0,0)*(- val(u.y,0,0,-1) +
         (val(u.z,0,1,-1) + val(u.z,0,1,0))/4. -
         (val(u.z,0,-1,-1) + val(u.z,0,-1,0))/4.)


      + val(mu.x,1,0,0)*(val(u.y,1,0,0) +
       (val(u.x,0,1,0) + val(u.x,1,1,0))/4. -
       (val(u.x,0,-1,0) + val(u.x,1,-1,0))/4.)
      - val(mu.x,0,0,0)*(- val(u.y,-1,0,0) +
         (val(u.x,-1,1,0) + val(u.x,0,1,0))/4. -
         (val(u.x,-1,-1,0) + val(u.x,0,-1,0))/4.)

      ) + val(r.y,0,0,0)*sq(Delta))/
    (sq(Delta)*((coord){1.,1.,1.}).y + dt/val(rho,0,0,0)*(2.*val(mu.y,0,1,0) + 2.*val(mu.y,0,0,0)

          + val(mu.z,0,0,1) + val(mu.z,0,0,0)


          + val(mu.x,1,0,0) + val(mu.x,0,0,0)

        ));
      
#line 151
val(w.z,0,0,0) = (dt/val(rho,0,0,0)*(2.*val(mu.z,0,0,1)*val(u.z,0,0,1) + 2.*val(mu.z,0,0,0)*val(u.z,0,0,-1)

      + val(mu.x,1,0,0)*(val(u.z,1,0,0) +
     (val(u.x,0,0,1) + val(u.x,1,0,1))/4. -
     (val(u.x,0,0,-1) + val(u.x,1,0,-1))/4.)
      - val(mu.x,0,0,0)*(- val(u.z,-1,0,0) +
         (val(u.x,-1,0,1) + val(u.x,0,0,1))/4. -
         (val(u.x,-1,0,-1) + val(u.x,0,0,-1))/4.)


      + val(mu.y,0,1,0)*(val(u.z,0,1,0) +
       (val(u.y,0,0,1) + val(u.y,0,1,1))/4. -
       (val(u.y,0,0,-1) + val(u.y,0,1,-1))/4.)
      - val(mu.y,0,0,0)*(- val(u.z,0,-1,0) +
         (val(u.y,0,-1,1) + val(u.y,0,0,1))/4. -
         (val(u.y,0,-1,-1) + val(u.y,0,0,-1))/4.)

      ) + val(r.z,0,0,0)*sq(Delta))/
    (sq(Delta)*((coord){1.,1.,1.}).z + dt/val(rho,0,0,0)*(2.*val(mu.z,0,0,1) + 2.*val(mu.z,0,0,0)

          + val(mu.x,1,0,0) + val(mu.x,0,0,0)


          + val(mu.y,0,1,0) + val(mu.y,0,0,0)

        ));
  }end_foreach_level_or_leaf();}}else if(is_constant(rho) && !is_constant(mu.x)){double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);


  {
#line 149
foreach_level_or_leaf (l) {
    
      val(w.x,0,0,0) = (dt/_const_rho*(2.*val(mu.x,1,0,0)*val(u.x,1,0,0) + 2.*val(mu.x,0,0,0)*val(u.x,-1,0,0)

      + val(mu.y,0,1,0)*(val(u.x,0,1,0) +
     (val(u.y,1,0,0) + val(u.y,1,1,0))/4. -
     (val(u.y,-1,0,0) + val(u.y,-1,1,0))/4.)
      - val(mu.y,0,0,0)*(- val(u.x,0,-1,0) +
         (val(u.y,1,-1,0) + val(u.y,1,0,0))/4. -
         (val(u.y,-1,-1,0) + val(u.y,-1,0,0))/4.)


      + val(mu.z,0,0,1)*(val(u.x,0,0,1) +
       (val(u.z,1,0,0) + val(u.z,1,0,1))/4. -
       (val(u.z,-1,0,0) + val(u.z,-1,0,1))/4.)
      - val(mu.z,0,0,0)*(- val(u.x,0,0,-1) +
         (val(u.z,1,0,-1) + val(u.z,1,0,0))/4. -
         (val(u.z,-1,0,-1) + val(u.z,-1,0,0))/4.)

      ) + val(r.x,0,0,0)*sq(Delta))/
    (sq(Delta)*((coord){1.,1.,1.}).x + dt/_const_rho*(2.*val(mu.x,1,0,0) + 2.*val(mu.x,0,0,0)

          + val(mu.y,0,1,0) + val(mu.y,0,0,0)


          + val(mu.z,0,0,1) + val(mu.z,0,0,0)

        ));
      
#line 151
val(w.y,0,0,0) = (dt/_const_rho*(2.*val(mu.y,0,1,0)*val(u.y,0,1,0) + 2.*val(mu.y,0,0,0)*val(u.y,0,-1,0)

      + val(mu.z,0,0,1)*(val(u.y,0,0,1) +
     (val(u.z,0,1,0) + val(u.z,0,1,1))/4. -
     (val(u.z,0,-1,0) + val(u.z,0,-1,1))/4.)
      - val(mu.z,0,0,0)*(- val(u.y,0,0,-1) +
         (val(u.z,0,1,-1) + val(u.z,0,1,0))/4. -
         (val(u.z,0,-1,-1) + val(u.z,0,-1,0))/4.)


      + val(mu.x,1,0,0)*(val(u.y,1,0,0) +
       (val(u.x,0,1,0) + val(u.x,1,1,0))/4. -
       (val(u.x,0,-1,0) + val(u.x,1,-1,0))/4.)
      - val(mu.x,0,0,0)*(- val(u.y,-1,0,0) +
         (val(u.x,-1,1,0) + val(u.x,0,1,0))/4. -
         (val(u.x,-1,-1,0) + val(u.x,0,-1,0))/4.)

      ) + val(r.y,0,0,0)*sq(Delta))/
    (sq(Delta)*((coord){1.,1.,1.}).y + dt/_const_rho*(2.*val(mu.y,0,1,0) + 2.*val(mu.y,0,0,0)

          + val(mu.z,0,0,1) + val(mu.z,0,0,0)


          + val(mu.x,1,0,0) + val(mu.x,0,0,0)

        ));
      
#line 151
val(w.z,0,0,0) = (dt/_const_rho*(2.*val(mu.z,0,0,1)*val(u.z,0,0,1) + 2.*val(mu.z,0,0,0)*val(u.z,0,0,-1)

      + val(mu.x,1,0,0)*(val(u.z,1,0,0) +
     (val(u.x,0,0,1) + val(u.x,1,0,1))/4. -
     (val(u.x,0,0,-1) + val(u.x,1,0,-1))/4.)
      - val(mu.x,0,0,0)*(- val(u.z,-1,0,0) +
         (val(u.x,-1,0,1) + val(u.x,0,0,1))/4. -
         (val(u.x,-1,0,-1) + val(u.x,0,0,-1))/4.)


      + val(mu.y,0,1,0)*(val(u.z,0,1,0) +
       (val(u.y,0,0,1) + val(u.y,0,1,1))/4. -
       (val(u.y,0,0,-1) + val(u.y,0,1,-1))/4.)
      - val(mu.y,0,0,0)*(- val(u.z,0,-1,0) +
         (val(u.y,0,-1,1) + val(u.y,0,0,1))/4. -
         (val(u.y,0,-1,-1) + val(u.y,0,0,-1))/4.)

      ) + val(r.z,0,0,0)*sq(Delta))/
    (sq(Delta)*((coord){1.,1.,1.}).z + dt/_const_rho*(2.*val(mu.z,0,0,1) + 2.*val(mu.z,0,0,0)

          + val(mu.x,1,0,0) + val(mu.x,0,0,0)


          + val(mu.y,0,1,0) + val(mu.y,0,0,0)

        ));
  }end_foreach_level_or_leaf();}}else if(!is_constant(rho) && is_constant(mu.x)){_coord _const_mu={_constant[mu.x.i-_NVARMAX],_constant[mu.y.i-_NVARMAX],_constant[mu.z.i-_NVARMAX]};NOT_UNUSED(_const_mu);


  {
#line 149
foreach_level_or_leaf (l) {
    
      val(w.x,0,0,0) = (dt/val(rho,0,0,0)*(2.*_const_mu.x*val(u.x,1,0,0) + 2.*_const_mu.x*val(u.x,-1,0,0)

      + _const_mu.y*(val(u.x,0,1,0) +
     (val(u.y,1,0,0) + val(u.y,1,1,0))/4. -
     (val(u.y,-1,0,0) + val(u.y,-1,1,0))/4.)
      - _const_mu.y*(- val(u.x,0,-1,0) +
         (val(u.y,1,-1,0) + val(u.y,1,0,0))/4. -
         (val(u.y,-1,-1,0) + val(u.y,-1,0,0))/4.)


      + _const_mu.z*(val(u.x,0,0,1) +
       (val(u.z,1,0,0) + val(u.z,1,0,1))/4. -
       (val(u.z,-1,0,0) + val(u.z,-1,0,1))/4.)
      - _const_mu.z*(- val(u.x,0,0,-1) +
         (val(u.z,1,0,-1) + val(u.z,1,0,0))/4. -
         (val(u.z,-1,0,-1) + val(u.z,-1,0,0))/4.)

      ) + val(r.x,0,0,0)*sq(Delta))/
    (sq(Delta)*((coord){1.,1.,1.}).x + dt/val(rho,0,0,0)*(2.*_const_mu.x + 2.*_const_mu.x

          + _const_mu.y + _const_mu.y


          + _const_mu.z + _const_mu.z

        ));
      
#line 151
val(w.y,0,0,0) = (dt/val(rho,0,0,0)*(2.*_const_mu.y*val(u.y,0,1,0) + 2.*_const_mu.y*val(u.y,0,-1,0)

      + _const_mu.z*(val(u.y,0,0,1) +
     (val(u.z,0,1,0) + val(u.z,0,1,1))/4. -
     (val(u.z,0,-1,0) + val(u.z,0,-1,1))/4.)
      - _const_mu.z*(- val(u.y,0,0,-1) +
         (val(u.z,0,1,-1) + val(u.z,0,1,0))/4. -
         (val(u.z,0,-1,-1) + val(u.z,0,-1,0))/4.)


      + _const_mu.x*(val(u.y,1,0,0) +
       (val(u.x,0,1,0) + val(u.x,1,1,0))/4. -
       (val(u.x,0,-1,0) + val(u.x,1,-1,0))/4.)
      - _const_mu.x*(- val(u.y,-1,0,0) +
         (val(u.x,-1,1,0) + val(u.x,0,1,0))/4. -
         (val(u.x,-1,-1,0) + val(u.x,0,-1,0))/4.)

      ) + val(r.y,0,0,0)*sq(Delta))/
    (sq(Delta)*((coord){1.,1.,1.}).y + dt/val(rho,0,0,0)*(2.*_const_mu.y + 2.*_const_mu.y

          + _const_mu.z + _const_mu.z


          + _const_mu.x + _const_mu.x

        ));
      
#line 151
val(w.z,0,0,0) = (dt/val(rho,0,0,0)*(2.*_const_mu.z*val(u.z,0,0,1) + 2.*_const_mu.z*val(u.z,0,0,-1)

      + _const_mu.x*(val(u.z,1,0,0) +
     (val(u.x,0,0,1) + val(u.x,1,0,1))/4. -
     (val(u.x,0,0,-1) + val(u.x,1,0,-1))/4.)
      - _const_mu.x*(- val(u.z,-1,0,0) +
         (val(u.x,-1,0,1) + val(u.x,0,0,1))/4. -
         (val(u.x,-1,0,-1) + val(u.x,0,0,-1))/4.)


      + _const_mu.y*(val(u.z,0,1,0) +
       (val(u.y,0,0,1) + val(u.y,0,1,1))/4. -
       (val(u.y,0,0,-1) + val(u.y,0,1,-1))/4.)
      - _const_mu.y*(- val(u.z,0,-1,0) +
         (val(u.y,0,-1,1) + val(u.y,0,0,1))/4. -
         (val(u.y,0,-1,-1) + val(u.y,0,0,-1))/4.)

      ) + val(r.z,0,0,0)*sq(Delta))/
    (sq(Delta)*((coord){1.,1.,1.}).z + dt/val(rho,0,0,0)*(2.*_const_mu.z + 2.*_const_mu.z

          + _const_mu.x + _const_mu.x


          + _const_mu.y + _const_mu.y

        ));
  }end_foreach_level_or_leaf();}}else {double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);_coord _const_mu={_constant[mu.x.i-_NVARMAX],_constant[mu.y.i-_NVARMAX],_constant[mu.z.i-_NVARMAX]};NOT_UNUSED(_const_mu);


  {
#line 149
foreach_level_or_leaf (l) {
    
      val(w.x,0,0,0) = (dt/_const_rho*(2.*_const_mu.x*val(u.x,1,0,0) + 2.*_const_mu.x*val(u.x,-1,0,0)

      + _const_mu.y*(val(u.x,0,1,0) +
     (val(u.y,1,0,0) + val(u.y,1,1,0))/4. -
     (val(u.y,-1,0,0) + val(u.y,-1,1,0))/4.)
      - _const_mu.y*(- val(u.x,0,-1,0) +
         (val(u.y,1,-1,0) + val(u.y,1,0,0))/4. -
         (val(u.y,-1,-1,0) + val(u.y,-1,0,0))/4.)


      + _const_mu.z*(val(u.x,0,0,1) +
       (val(u.z,1,0,0) + val(u.z,1,0,1))/4. -
       (val(u.z,-1,0,0) + val(u.z,-1,0,1))/4.)
      - _const_mu.z*(- val(u.x,0,0,-1) +
         (val(u.z,1,0,-1) + val(u.z,1,0,0))/4. -
         (val(u.z,-1,0,-1) + val(u.z,-1,0,0))/4.)

      ) + val(r.x,0,0,0)*sq(Delta))/
    (sq(Delta)*((coord){1.,1.,1.}).x + dt/_const_rho*(2.*_const_mu.x + 2.*_const_mu.x

          + _const_mu.y + _const_mu.y


          + _const_mu.z + _const_mu.z

        ));
      
#line 151
val(w.y,0,0,0) = (dt/_const_rho*(2.*_const_mu.y*val(u.y,0,1,0) + 2.*_const_mu.y*val(u.y,0,-1,0)

      + _const_mu.z*(val(u.y,0,0,1) +
     (val(u.z,0,1,0) + val(u.z,0,1,1))/4. -
     (val(u.z,0,-1,0) + val(u.z,0,-1,1))/4.)
      - _const_mu.z*(- val(u.y,0,0,-1) +
         (val(u.z,0,1,-1) + val(u.z,0,1,0))/4. -
         (val(u.z,0,-1,-1) + val(u.z,0,-1,0))/4.)


      + _const_mu.x*(val(u.y,1,0,0) +
       (val(u.x,0,1,0) + val(u.x,1,1,0))/4. -
       (val(u.x,0,-1,0) + val(u.x,1,-1,0))/4.)
      - _const_mu.x*(- val(u.y,-1,0,0) +
         (val(u.x,-1,1,0) + val(u.x,0,1,0))/4. -
         (val(u.x,-1,-1,0) + val(u.x,0,-1,0))/4.)

      ) + val(r.y,0,0,0)*sq(Delta))/
    (sq(Delta)*((coord){1.,1.,1.}).y + dt/_const_rho*(2.*_const_mu.y + 2.*_const_mu.y

          + _const_mu.z + _const_mu.z


          + _const_mu.x + _const_mu.x

        ));
      
#line 151
val(w.z,0,0,0) = (dt/_const_rho*(2.*_const_mu.z*val(u.z,0,0,1) + 2.*_const_mu.z*val(u.z,0,0,-1)

      + _const_mu.x*(val(u.z,1,0,0) +
     (val(u.x,0,0,1) + val(u.x,1,0,1))/4. -
     (val(u.x,0,0,-1) + val(u.x,1,0,-1))/4.)
      - _const_mu.x*(- val(u.z,-1,0,0) +
         (val(u.x,-1,0,1) + val(u.x,0,0,1))/4. -
         (val(u.x,-1,0,-1) + val(u.x,0,0,-1))/4.)


      + _const_mu.y*(val(u.z,0,1,0) +
       (val(u.y,0,0,1) + val(u.y,0,1,1))/4. -
       (val(u.y,0,0,-1) + val(u.y,0,1,-1))/4.)
      - _const_mu.y*(- val(u.z,0,-1,0) +
         (val(u.y,0,-1,1) + val(u.y,0,0,1))/4. -
         (val(u.y,0,-1,-1) + val(u.y,0,0,-1))/4.)

      ) + val(r.z,0,0,0)*sq(Delta))/
    (sq(Delta)*((coord){1.,1.,1.}).z + dt/_const_rho*(2.*_const_mu.z + 2.*_const_mu.z

          + _const_mu.x + _const_mu.x


          + _const_mu.y + _const_mu.y

        ));
  }end_foreach_level_or_leaf();}}
#line 195 "/home/lisergey/basilisk/src/viscosity.h"
}
#line 204 "/home/lisergey/basilisk/src/viscosity.h"
static double residual_viscosity (scalar * a, scalar * b, scalar * resl,
      void * data)
{
  struct Viscosity * p = (struct Viscosity *) data;
          vector mu = p->mu;
          scalar rho = p->rho;
  double dt = p->dt;
  vector u = (*((vector *)&(a[0]))), r = (*((vector *)&(b[0]))), res = (*((vector *)&(resl[0])));
  double maxres = 0.;
#line 221 "/home/lisergey/basilisk/src/viscosity.h"
  boundary_internal ((scalar *)((vector[]){u,{{-1},{-1},{-1}}}), "/home/lisergey/basilisk/src/viscosity.h", 221);

   {
    vector  taux=new_face_vector("taux");
    if(!is_constant(mu.x)){
    
#line 225
foreach_face_stencil(1,{(NonLocal[]){{"u","vector",(void *)&u,NULL,0},{"mu","vector",(void *)&mu,NULL,0},{"taux","vector",(void *)&taux,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n// #line 225 \"/home/lisergey/basilisk/src/viscosity.h\"\nis_face_x(){\n      val_out_(taux.x,0,0,0) = 2.*val(mu.x,0,0,0)*(val(u.x,0,0,0) - val(u.x,-1,0,0))/Delta;}end_is_face_x()"})_stencil_is_face_x(){
      {_stencil_val(mu.x,0,0,0);_stencil_val(u.x,0,0,0); _stencil_val(u.x,-1,0,0);_stencil_val_a(taux.x,0,0,0);   }}end__stencil_is_face_x()end_foreach_face_stencil();{
#line 225
foreach_face_generic()is_face_x(){
      val(taux.x,0,0,0) = 2.*val(mu.x,0,0,0)*(val(u.x,0,0,0) - val(u.x,-1,0,0))/Delta;}end_is_face_x()end_foreach_face_generic();}}else {_coord _const_mu={_constant[mu.x.i-_NVARMAX],_constant[mu.y.i-_NVARMAX],_constant[mu.z.i-_NVARMAX]};NOT_UNUSED(_const_mu);
    
#line 225
foreach_face_stencil(1,{(NonLocal[]){{"u","vector",(void *)&u,NULL,0},{"_const_mu","_coord",(void *)&_const_mu,NULL,0},{"taux","vector",(void *)&taux,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n// #line 225 \"/home/lisergey/basilisk/src/viscosity.h\"\nis_face_x(){\n      val_out_(taux.x,0,0,0) = 2.*_const_mu.x*(val(u.x,0,0,0) - val(u.x,-1,0,0))/Delta;}end_is_face_x()"})_stencil_is_face_x(){
      {;_stencil_val(u.x,0,0,0); _stencil_val(u.x,-1,0,0);_stencil_val_a(taux.x,0,0,0);   }}end__stencil_is_face_x()end_foreach_face_stencil();
    {
#line 225
foreach_face_generic()is_face_x(){
      val(taux.x,0,0,0) = 2.*_const_mu.x*(val(u.x,0,0,0) - val(u.x,-1,0,0))/Delta;}end_is_face_x()end_foreach_face_generic();}}

      if(!is_constant(mu.x)){

      
#line 228
foreach_face_stencil(1,{(NonLocal[]){{"u","vector",(void *)&u,NULL,0},{"mu","vector",(void *)&mu,NULL,0},{"taux","vector",(void *)&taux,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n// #line 228 \"/home/lisergey/basilisk/src/viscosity.h\"\nis_face_y(){\n val_out_(taux.y,0,0,0) = val(mu.y,0,0,0)*(val(u.x,0,0,0) - val(u.x,0,-1,0) +\n      (val(u.y,1,-1,0) + val(u.y,1,0,0))/4. -\n      (val(u.y,-1,-1,0) + val(u.y,-1,0,0))/4.)/Delta;}end_is_face_y()"})_stencil_is_face_y(){
 { _stencil_val(mu.y,0,0,0);_stencil_val(u.x,0,0,0); _stencil_val(u.x,0,-1,0);
_stencil_val(u.y,1,-1,0); _stencil_val(u.y,1,0,0);
_stencil_val(u.y,-1,-1,0); _stencil_val(u.y,-1,0,0);
#line 229
_stencil_val_a(taux.y,0,0,0);    
        
      }}end__stencil_is_face_y()end_foreach_face_stencil();{
#line 228
foreach_face_generic()is_face_y(){
 val(taux.y,0,0,0) = val(mu.y,0,0,0)*(val(u.x,0,0,0) - val(u.x,0,-1,0) +
      (val(u.y,1,-1,0) + val(u.y,1,0,0))/4. -
      (val(u.y,-1,-1,0) + val(u.y,-1,0,0))/4.)/Delta;}end_is_face_y()end_foreach_face_generic();}}else {_coord _const_mu={_constant[mu.x.i-_NVARMAX],_constant[mu.y.i-_NVARMAX],_constant[mu.z.i-_NVARMAX]};NOT_UNUSED(_const_mu);

      
#line 228
foreach_face_stencil(1,{(NonLocal[]){{"u","vector",(void *)&u,NULL,0},{"_const_mu","_coord",(void *)&_const_mu,NULL,0},{"taux","vector",(void *)&taux,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n// #line 228 \"/home/lisergey/basilisk/src/viscosity.h\"\nis_face_y(){\n val_out_(taux.y,0,0,0) = _const_mu.y*(val(u.x,0,0,0) - val(u.x,0,-1,0) +\n      (val(u.y,1,-1,0) + val(u.y,1,0,0))/4. -\n      (val(u.y,-1,-1,0) + val(u.y,-1,0,0))/4.)/Delta;}end_is_face_y()"})_stencil_is_face_y(){
 {;_stencil_val(u.x,0,0,0); _stencil_val(u.x,0,-1,0);
_stencil_val(u.y,1,-1,0); _stencil_val(u.y,1,0,0);
_stencil_val(u.y,-1,-1,0); _stencil_val(u.y,-1,0,0);
#line 229
_stencil_val_a(taux.y,0,0,0);    
        
      }}end__stencil_is_face_y()end_foreach_face_stencil();

      {
#line 228
foreach_face_generic()is_face_y(){
 val(taux.y,0,0,0) = _const_mu.y*(val(u.x,0,0,0) - val(u.x,0,-1,0) +
      (val(u.y,1,-1,0) + val(u.y,1,0,0))/4. -
      (val(u.y,-1,-1,0) + val(u.y,-1,0,0))/4.)/Delta;}end_is_face_y()end_foreach_face_generic();}}


      if(!is_constant(mu.x)){


      
#line 234
foreach_face_stencil(1,{(NonLocal[]){{"u","vector",(void *)&u,NULL,0},{"mu","vector",(void *)&mu,NULL,0},{"taux","vector",(void *)&taux,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n// #line 234 \"/home/lisergey/basilisk/src/viscosity.h\"\nis_face_z(){\n val_out_(taux.z,0,0,0) = val(mu.z,0,0,0)*(val(u.x,0,0,0) - val(u.x,0,0,-1) +\n      (val(u.z,1,0,-1) + val(u.z,1,0,0))/4. -\n      (val(u.z,-1,0,-1) + val(u.z,-1,0,0))/4.)/Delta;}end_is_face_z()"})_stencil_is_face_z(){
 { _stencil_val(mu.z,0,0,0);_stencil_val(u.x,0,0,0); _stencil_val(u.x,0,0,-1);
_stencil_val(u.z,1,0,-1); _stencil_val(u.z,1,0,0);
_stencil_val(u.z,-1,0,-1); _stencil_val(u.z,-1,0,0);
#line 235
_stencil_val_a(taux.z,0,0,0);    
        
      }}end__stencil_is_face_z()end_foreach_face_stencil();{
#line 234
foreach_face_generic()is_face_z(){
 val(taux.z,0,0,0) = val(mu.z,0,0,0)*(val(u.x,0,0,0) - val(u.x,0,0,-1) +
      (val(u.z,1,0,-1) + val(u.z,1,0,0))/4. -
      (val(u.z,-1,0,-1) + val(u.z,-1,0,0))/4.)/Delta;}end_is_face_z()end_foreach_face_generic();}}else {_coord _const_mu={_constant[mu.x.i-_NVARMAX],_constant[mu.y.i-_NVARMAX],_constant[mu.z.i-_NVARMAX]};NOT_UNUSED(_const_mu);


      
#line 234
foreach_face_stencil(1,{(NonLocal[]){{"u","vector",(void *)&u,NULL,0},{"_const_mu","_coord",(void *)&_const_mu,NULL,0},{"taux","vector",(void *)&taux,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n// #line 234 \"/home/lisergey/basilisk/src/viscosity.h\"\nis_face_z(){\n val_out_(taux.z,0,0,0) = _const_mu.z*(val(u.x,0,0,0) - val(u.x,0,0,-1) +\n      (val(u.z,1,0,-1) + val(u.z,1,0,0))/4. -\n      (val(u.z,-1,0,-1) + val(u.z,-1,0,0))/4.)/Delta;}end_is_face_z()"})_stencil_is_face_z(){
 {;_stencil_val(u.x,0,0,0); _stencil_val(u.x,0,0,-1);
_stencil_val(u.z,1,0,-1); _stencil_val(u.z,1,0,0);
_stencil_val(u.z,-1,0,-1); _stencil_val(u.z,-1,0,0);
#line 235
_stencil_val_a(taux.z,0,0,0);    
        
      }}end__stencil_is_face_z()end_foreach_face_stencil();


      {
#line 234
foreach_face_generic()is_face_z(){
 val(taux.z,0,0,0) = _const_mu.z*(val(u.x,0,0,0) - val(u.x,0,0,-1) +
      (val(u.z,1,0,-1) + val(u.z,1,0,0))/4. -
      (val(u.z,-1,0,-1) + val(u.z,-1,0,0))/4.)/Delta;}end_is_face_z()end_foreach_face_generic();}}

    if(!is_constant(rho)){

    
#line 239
foreach_stencil (1,{(NonLocal[]){{"maxres","double",(void *)&maxres,NULL,0,'M'},{"rho","scalar",(void *)&rho,NULL,0},{"dt","double",(void *)&dt,NULL,0},{"u","vector",(void *)&u,NULL,0},{"r","vector",(void *)&r,NULL,0},{"res","vector",(void *)&res,NULL,0},{"taux","vector",(void *)&taux,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},""," \n// #line 239 \"/home/lisergey/basilisk/src/viscosity.h\"\n{\n      real d = 0.;\n      \n d += val(taux.x,1,0,0) - val(taux.x,0,0,0);\n \n// #line 242\nd += val(taux.y,0,1,0) - val(taux.y,0,0,0);\n \n// #line 242\nd += val(taux.z,0,0,1) - val(taux.z,0,0,0);\n      val_out_(res.x,0,0,0) = val(r.x,0,0,0) - ((coord){1.,1.,1.}).x*val(u.x,0,0,0) + dt/val(rho,0,0,0)*d/Delta;\n      if (fabs (val(res.x,0,0,0)) > maxres)\n maxres = fabs (val(res.x,0,0,0));\n    }"}) {   
      
      
 { _stencil_val(taux.x,1,0,0); _stencil_val(taux.x,0,0,0);  }
 
#line 242
{ _stencil_val(taux.y,0,1,0); _stencil_val(taux.y,0,0,0);  }
 
#line 242
{ _stencil_val(taux.z,0,0,1); _stencil_val(taux.z,0,0,0);  } 
_stencil_val(r.x,0,0,0);_stencil_val(u.x,0,0,0);_stencil_val(rho,0,0,0);
      
#line 243
_stencil_val_a(res.x,0,0,0);
_stencil_val(res.x,0,0,0);
 {_stencil_val(res.x,0,0,0);   }     
          
    
#line 246
}end_foreach_stencil();
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel reduction(max:maxres)){
#line 239
foreach () {
      double d = 0.;
      
 d += val(taux.x,1,0,0) - val(taux.x,0,0,0);
 
#line 242
d += val(taux.y,0,1,0) - val(taux.y,0,0,0);
 
#line 242
d += val(taux.z,0,0,1) - val(taux.z,0,0,0);
      val(res.x,0,0,0) = val(r.x,0,0,0) - ((coord){1.,1.,1.}).x*val(u.x,0,0,0) + dt/val(rho,0,0,0)*d/Delta;
      if (fabs (val(res.x,0,0,0)) > maxres)
 maxres = fabs (val(res.x,0,0,0));
    }end_foreach();mpi_all_reduce_array(&maxres,double,MPI_MAX,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 246
}else {double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);

    
#line 239
foreach_stencil (1,{(NonLocal[]){{"maxres","double",(void *)&maxres,NULL,0,'M'},{"_const_rho","double",(void *)&_const_rho,NULL,0},{"dt","double",(void *)&dt,NULL,0},{"u","vector",(void *)&u,NULL,0},{"r","vector",(void *)&r,NULL,0},{"res","vector",(void *)&res,NULL,0},{"taux","vector",(void *)&taux,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},""," \n// #line 239 \"/home/lisergey/basilisk/src/viscosity.h\"\n{\n      real d = 0.;\n      \n d += val(taux.x,1,0,0) - val(taux.x,0,0,0);\n \n// #line 242\nd += val(taux.y,0,1,0) - val(taux.y,0,0,0);\n \n// #line 242\nd += val(taux.z,0,0,1) - val(taux.z,0,0,0);\n      val_out_(res.x,0,0,0) = val(r.x,0,0,0) - ((coord){1.,1.,1.}).x*val(u.x,0,0,0) + dt/_const_rho*d/Delta;\n      if (fabs (val(res.x,0,0,0)) > maxres)\n maxres = fabs (val(res.x,0,0,0));\n    }"}) {   
      
      
 { _stencil_val(taux.x,1,0,0); _stencil_val(taux.x,0,0,0);  }
 
#line 242
{ _stencil_val(taux.y,0,1,0); _stencil_val(taux.y,0,0,0);  }
 
#line 242
{ _stencil_val(taux.z,0,0,1); _stencil_val(taux.z,0,0,0);  } 
_stencil_val(r.x,0,0,0);_stencil_val(u.x,0,0,0);;
      
#line 243
_stencil_val_a(res.x,0,0,0);
_stencil_val(res.x,0,0,0);
 {_stencil_val(res.x,0,0,0);   }     
          
    
#line 246
}end_foreach_stencil();

    
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel reduction(max:maxres)){
#line 239
foreach () {
      double d = 0.;
      
 d += val(taux.x,1,0,0) - val(taux.x,0,0,0);
 
#line 242
d += val(taux.y,0,1,0) - val(taux.y,0,0,0);
 
#line 242
d += val(taux.z,0,0,1) - val(taux.z,0,0,0);
      val(res.x,0,0,0) = val(r.x,0,0,0) - ((coord){1.,1.,1.}).x*val(u.x,0,0,0) + dt/_const_rho*d/Delta;
      if (fabs (val(res.x,0,0,0)) > maxres)
 maxres = fabs (val(res.x,0,0,0));
    }end_foreach();mpi_all_reduce_array(&maxres,double,MPI_MAX,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 246
}delete((scalar*)((vector[]){taux,{{-1},{-1},{-1}}}));
  } 
#line 223
{
    vector  taux=new_face_vector("taux");
    if(!is_constant(mu.y)){
    
#line 225
foreach_face_stencil(1,{(NonLocal[]){{"u","vector",(void *)&u,NULL,0},{"mu","vector",(void *)&mu,NULL,0},{"taux","vector",(void *)&taux,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n// #line 225 \"/home/lisergey/basilisk/src/viscosity.h\"\nis_face_y(){\n      val_out_(taux.y,0,0,0) = 2.*val(mu.y,0,0,0)*(val(u.y,0,0,0) - val(u.y,0,-1,0))/Delta;}end_is_face_y()"})_stencil_is_face_y(){
      {_stencil_val(mu.y,0,0,0);_stencil_val(u.y,0,0,0); _stencil_val(u.y,0,-1,0);_stencil_val_a(taux.y,0,0,0);   }}end__stencil_is_face_y()end_foreach_face_stencil();{
#line 225
foreach_face_generic()is_face_y(){
      val(taux.y,0,0,0) = 2.*val(mu.y,0,0,0)*(val(u.y,0,0,0) - val(u.y,0,-1,0))/Delta;}end_is_face_y()end_foreach_face_generic();}}else {_coord _const_mu={_constant[mu.y.i-_NVARMAX],_constant[mu.z.i-_NVARMAX],_constant[mu.x.i-_NVARMAX]};NOT_UNUSED(_const_mu);
    
#line 225
foreach_face_stencil(1,{(NonLocal[]){{"u","vector",(void *)&u,NULL,0},{"_const_mu","_coord",(void *)&_const_mu,NULL,0},{"taux","vector",(void *)&taux,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n// #line 225 \"/home/lisergey/basilisk/src/viscosity.h\"\nis_face_y(){\n      val_out_(taux.y,0,0,0) = 2.*_const_mu.y*(val(u.y,0,0,0) - val(u.y,0,-1,0))/Delta;}end_is_face_y()"})_stencil_is_face_y(){
      {;_stencil_val(u.y,0,0,0); _stencil_val(u.y,0,-1,0);_stencil_val_a(taux.y,0,0,0);   }}end__stencil_is_face_y()end_foreach_face_stencil();
    {
#line 225
foreach_face_generic()is_face_y(){
      val(taux.y,0,0,0) = 2.*_const_mu.y*(val(u.y,0,0,0) - val(u.y,0,-1,0))/Delta;}end_is_face_y()end_foreach_face_generic();}}

      if(!is_constant(mu.y)){

      
#line 228
foreach_face_stencil(1,{(NonLocal[]){{"u","vector",(void *)&u,NULL,0},{"mu","vector",(void *)&mu,NULL,0},{"taux","vector",(void *)&taux,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n// #line 228 \"/home/lisergey/basilisk/src/viscosity.h\"\nis_face_z(){\n val_out_(taux.z,0,0,0) = val(mu.z,0,0,0)*(val(u.y,0,0,0) - val(u.y,0,0,-1) +\n      (val(u.z,0,1,-1) + val(u.z,0,1,0))/4. -\n      (val(u.z,0,-1,-1) + val(u.z,0,-1,0))/4.)/Delta;}end_is_face_z()"})_stencil_is_face_z(){
 { _stencil_val(mu.z,0,0,0);_stencil_val(u.y,0,0,0); _stencil_val(u.y,0,0,-1);
_stencil_val(u.z,0,1,-1); _stencil_val(u.z,0,1,0);
_stencil_val(u.z,0,-1,-1); _stencil_val(u.z,0,-1,0);
#line 229
_stencil_val_a(taux.z,0,0,0);    
        
      }}end__stencil_is_face_z()end_foreach_face_stencil();{
#line 228
foreach_face_generic()is_face_z(){
 val(taux.z,0,0,0) = val(mu.z,0,0,0)*(val(u.y,0,0,0) - val(u.y,0,0,-1) +
      (val(u.z,0,1,-1) + val(u.z,0,1,0))/4. -
      (val(u.z,0,-1,-1) + val(u.z,0,-1,0))/4.)/Delta;}end_is_face_z()end_foreach_face_generic();}}else {_coord _const_mu={_constant[mu.y.i-_NVARMAX],_constant[mu.z.i-_NVARMAX],_constant[mu.x.i-_NVARMAX]};NOT_UNUSED(_const_mu);

      
#line 228
foreach_face_stencil(1,{(NonLocal[]){{"u","vector",(void *)&u,NULL,0},{"_const_mu","_coord",(void *)&_const_mu,NULL,0},{"taux","vector",(void *)&taux,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n// #line 228 \"/home/lisergey/basilisk/src/viscosity.h\"\nis_face_z(){\n val_out_(taux.z,0,0,0) = _const_mu.z*(val(u.y,0,0,0) - val(u.y,0,0,-1) +\n      (val(u.z,0,1,-1) + val(u.z,0,1,0))/4. -\n      (val(u.z,0,-1,-1) + val(u.z,0,-1,0))/4.)/Delta;}end_is_face_z()"})_stencil_is_face_z(){
 {;_stencil_val(u.y,0,0,0); _stencil_val(u.y,0,0,-1);
_stencil_val(u.z,0,1,-1); _stencil_val(u.z,0,1,0);
_stencil_val(u.z,0,-1,-1); _stencil_val(u.z,0,-1,0);
#line 229
_stencil_val_a(taux.z,0,0,0);    
        
      }}end__stencil_is_face_z()end_foreach_face_stencil();

      {
#line 228
foreach_face_generic()is_face_z(){
 val(taux.z,0,0,0) = _const_mu.z*(val(u.y,0,0,0) - val(u.y,0,0,-1) +
      (val(u.z,0,1,-1) + val(u.z,0,1,0))/4. -
      (val(u.z,0,-1,-1) + val(u.z,0,-1,0))/4.)/Delta;}end_is_face_z()end_foreach_face_generic();}}


      if(!is_constant(mu.y)){


      
#line 234
foreach_face_stencil(1,{(NonLocal[]){{"u","vector",(void *)&u,NULL,0},{"mu","vector",(void *)&mu,NULL,0},{"taux","vector",(void *)&taux,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n// #line 234 \"/home/lisergey/basilisk/src/viscosity.h\"\nis_face_x(){\n val_out_(taux.x,0,0,0) = val(mu.x,0,0,0)*(val(u.y,0,0,0) - val(u.y,-1,0,0) +\n      (val(u.x,-1,1,0) + val(u.x,0,1,0))/4. -\n      (val(u.x,-1,-1,0) + val(u.x,0,-1,0))/4.)/Delta;}end_is_face_x()"})_stencil_is_face_x(){
 { _stencil_val(mu.x,0,0,0);_stencil_val(u.y,0,0,0); _stencil_val(u.y,-1,0,0);
_stencil_val(u.x,-1,1,0); _stencil_val(u.x,0,1,0);
_stencil_val(u.x,-1,-1,0); _stencil_val(u.x,0,-1,0);
#line 235
_stencil_val_a(taux.x,0,0,0);    
        
      }}end__stencil_is_face_x()end_foreach_face_stencil();{
#line 234
foreach_face_generic()is_face_x(){
 val(taux.x,0,0,0) = val(mu.x,0,0,0)*(val(u.y,0,0,0) - val(u.y,-1,0,0) +
      (val(u.x,-1,1,0) + val(u.x,0,1,0))/4. -
      (val(u.x,-1,-1,0) + val(u.x,0,-1,0))/4.)/Delta;}end_is_face_x()end_foreach_face_generic();}}else {_coord _const_mu={_constant[mu.y.i-_NVARMAX],_constant[mu.z.i-_NVARMAX],_constant[mu.x.i-_NVARMAX]};NOT_UNUSED(_const_mu);


      
#line 234
foreach_face_stencil(1,{(NonLocal[]){{"u","vector",(void *)&u,NULL,0},{"_const_mu","_coord",(void *)&_const_mu,NULL,0},{"taux","vector",(void *)&taux,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n// #line 234 \"/home/lisergey/basilisk/src/viscosity.h\"\nis_face_x(){\n val_out_(taux.x,0,0,0) = _const_mu.x*(val(u.y,0,0,0) - val(u.y,-1,0,0) +\n      (val(u.x,-1,1,0) + val(u.x,0,1,0))/4. -\n      (val(u.x,-1,-1,0) + val(u.x,0,-1,0))/4.)/Delta;}end_is_face_x()"})_stencil_is_face_x(){
 {;_stencil_val(u.y,0,0,0); _stencil_val(u.y,-1,0,0);
_stencil_val(u.x,-1,1,0); _stencil_val(u.x,0,1,0);
_stencil_val(u.x,-1,-1,0); _stencil_val(u.x,0,-1,0);
#line 235
_stencil_val_a(taux.x,0,0,0);    
        
      }}end__stencil_is_face_x()end_foreach_face_stencil();


      {
#line 234
foreach_face_generic()is_face_x(){
 val(taux.x,0,0,0) = _const_mu.x*(val(u.y,0,0,0) - val(u.y,-1,0,0) +
      (val(u.x,-1,1,0) + val(u.x,0,1,0))/4. -
      (val(u.x,-1,-1,0) + val(u.x,0,-1,0))/4.)/Delta;}end_is_face_x()end_foreach_face_generic();}}

    if(!is_constant(rho)){

    
#line 239
foreach_stencil (1,{(NonLocal[]){{"maxres","double",(void *)&maxres,NULL,0,'M'},{"rho","scalar",(void *)&rho,NULL,0},{"dt","double",(void *)&dt,NULL,0},{"u","vector",(void *)&u,NULL,0},{"r","vector",(void *)&r,NULL,0},{"res","vector",(void *)&res,NULL,0},{"taux","vector",(void *)&taux,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},""," \n// #line 239 \"/home/lisergey/basilisk/src/viscosity.h\"\n{\n      real d = 0.;\n      \n d += val(taux.y,0,1,0) - val(taux.y,0,0,0);\n \n// #line 242\nd += val(taux.z,0,0,1) - val(taux.z,0,0,0);\n \n// #line 242\nd += val(taux.x,1,0,0) - val(taux.x,0,0,0);\n      val_out_(res.y,0,0,0) = val(r.y,0,0,0) - ((coord){1.,1.,1.}).y*val(u.y,0,0,0) + dt/val(rho,0,0,0)*d/Delta;\n      if (fabs (val(res.y,0,0,0)) > maxres)\n maxres = fabs (val(res.y,0,0,0));\n    }"}) {   
      
      
 { _stencil_val(taux.y,0,1,0); _stencil_val(taux.y,0,0,0);  }
 
#line 242
{ _stencil_val(taux.z,0,0,1); _stencil_val(taux.z,0,0,0);  }
 
#line 242
{ _stencil_val(taux.x,1,0,0); _stencil_val(taux.x,0,0,0);  } 
_stencil_val(r.y,0,0,0);_stencil_val(u.y,0,0,0);_stencil_val(rho,0,0,0);
      
#line 243
_stencil_val_a(res.y,0,0,0);
_stencil_val(res.y,0,0,0);
 {_stencil_val(res.y,0,0,0);   }     
          
    
#line 246
}end_foreach_stencil();
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel reduction(max:maxres)){
#line 239
foreach () {
      double d = 0.;
      
 d += val(taux.y,0,1,0) - val(taux.y,0,0,0);
 
#line 242
d += val(taux.z,0,0,1) - val(taux.z,0,0,0);
 
#line 242
d += val(taux.x,1,0,0) - val(taux.x,0,0,0);
      val(res.y,0,0,0) = val(r.y,0,0,0) - ((coord){1.,1.,1.}).y*val(u.y,0,0,0) + dt/val(rho,0,0,0)*d/Delta;
      if (fabs (val(res.y,0,0,0)) > maxres)
 maxres = fabs (val(res.y,0,0,0));
    }end_foreach();mpi_all_reduce_array(&maxres,double,MPI_MAX,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 246
}else {double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);

    
#line 239
foreach_stencil (1,{(NonLocal[]){{"maxres","double",(void *)&maxres,NULL,0,'M'},{"_const_rho","double",(void *)&_const_rho,NULL,0},{"dt","double",(void *)&dt,NULL,0},{"u","vector",(void *)&u,NULL,0},{"r","vector",(void *)&r,NULL,0},{"res","vector",(void *)&res,NULL,0},{"taux","vector",(void *)&taux,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},""," \n// #line 239 \"/home/lisergey/basilisk/src/viscosity.h\"\n{\n      real d = 0.;\n      \n d += val(taux.y,0,1,0) - val(taux.y,0,0,0);\n \n// #line 242\nd += val(taux.z,0,0,1) - val(taux.z,0,0,0);\n \n// #line 242\nd += val(taux.x,1,0,0) - val(taux.x,0,0,0);\n      val_out_(res.y,0,0,0) = val(r.y,0,0,0) - ((coord){1.,1.,1.}).y*val(u.y,0,0,0) + dt/_const_rho*d/Delta;\n      if (fabs (val(res.y,0,0,0)) > maxres)\n maxres = fabs (val(res.y,0,0,0));\n    }"}) {   
      
      
 { _stencil_val(taux.y,0,1,0); _stencil_val(taux.y,0,0,0);  }
 
#line 242
{ _stencil_val(taux.z,0,0,1); _stencil_val(taux.z,0,0,0);  }
 
#line 242
{ _stencil_val(taux.x,1,0,0); _stencil_val(taux.x,0,0,0);  } 
_stencil_val(r.y,0,0,0);_stencil_val(u.y,0,0,0);;
      
#line 243
_stencil_val_a(res.y,0,0,0);
_stencil_val(res.y,0,0,0);
 {_stencil_val(res.y,0,0,0);   }     
          
    
#line 246
}end_foreach_stencil();

    
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel reduction(max:maxres)){
#line 239
foreach () {
      double d = 0.;
      
 d += val(taux.y,0,1,0) - val(taux.y,0,0,0);
 
#line 242
d += val(taux.z,0,0,1) - val(taux.z,0,0,0);
 
#line 242
d += val(taux.x,1,0,0) - val(taux.x,0,0,0);
      val(res.y,0,0,0) = val(r.y,0,0,0) - ((coord){1.,1.,1.}).y*val(u.y,0,0,0) + dt/_const_rho*d/Delta;
      if (fabs (val(res.y,0,0,0)) > maxres)
 maxres = fabs (val(res.y,0,0,0));
    }end_foreach();mpi_all_reduce_array(&maxres,double,MPI_MAX,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 246
}delete((scalar*)((vector[]){taux,{{-1},{-1},{-1}}}));
  } 
#line 223
{
    vector  taux=new_face_vector("taux");
    if(!is_constant(mu.z)){
    
#line 225
foreach_face_stencil(1,{(NonLocal[]){{"u","vector",(void *)&u,NULL,0},{"mu","vector",(void *)&mu,NULL,0},{"taux","vector",(void *)&taux,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n// #line 225 \"/home/lisergey/basilisk/src/viscosity.h\"\nis_face_z(){\n      val_out_(taux.z,0,0,0) = 2.*val(mu.z,0,0,0)*(val(u.z,0,0,0) - val(u.z,0,0,-1))/Delta;}end_is_face_z()"})_stencil_is_face_z(){
      {_stencil_val(mu.z,0,0,0);_stencil_val(u.z,0,0,0); _stencil_val(u.z,0,0,-1);_stencil_val_a(taux.z,0,0,0);   }}end__stencil_is_face_z()end_foreach_face_stencil();{
#line 225
foreach_face_generic()is_face_z(){
      val(taux.z,0,0,0) = 2.*val(mu.z,0,0,0)*(val(u.z,0,0,0) - val(u.z,0,0,-1))/Delta;}end_is_face_z()end_foreach_face_generic();}}else {_coord _const_mu={_constant[mu.z.i-_NVARMAX],_constant[mu.x.i-_NVARMAX],_constant[mu.y.i-_NVARMAX]};NOT_UNUSED(_const_mu);
    
#line 225
foreach_face_stencil(1,{(NonLocal[]){{"u","vector",(void *)&u,NULL,0},{"_const_mu","_coord",(void *)&_const_mu,NULL,0},{"taux","vector",(void *)&taux,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n// #line 225 \"/home/lisergey/basilisk/src/viscosity.h\"\nis_face_z(){\n      val_out_(taux.z,0,0,0) = 2.*_const_mu.z*(val(u.z,0,0,0) - val(u.z,0,0,-1))/Delta;}end_is_face_z()"})_stencil_is_face_z(){
      {;_stencil_val(u.z,0,0,0); _stencil_val(u.z,0,0,-1);_stencil_val_a(taux.z,0,0,0);   }}end__stencil_is_face_z()end_foreach_face_stencil();
    {
#line 225
foreach_face_generic()is_face_z(){
      val(taux.z,0,0,0) = 2.*_const_mu.z*(val(u.z,0,0,0) - val(u.z,0,0,-1))/Delta;}end_is_face_z()end_foreach_face_generic();}}

      if(!is_constant(mu.z)){

      
#line 228
foreach_face_stencil(1,{(NonLocal[]){{"u","vector",(void *)&u,NULL,0},{"mu","vector",(void *)&mu,NULL,0},{"taux","vector",(void *)&taux,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n// #line 228 \"/home/lisergey/basilisk/src/viscosity.h\"\nis_face_x(){\n val_out_(taux.x,0,0,0) = val(mu.x,0,0,0)*(val(u.z,0,0,0) - val(u.z,-1,0,0) +\n      (val(u.x,-1,0,1) + val(u.x,0,0,1))/4. -\n      (val(u.x,-1,0,-1) + val(u.x,0,0,-1))/4.)/Delta;}end_is_face_x()"})_stencil_is_face_x(){
 { _stencil_val(mu.x,0,0,0);_stencil_val(u.z,0,0,0); _stencil_val(u.z,-1,0,0);
_stencil_val(u.x,-1,0,1); _stencil_val(u.x,0,0,1);
_stencil_val(u.x,-1,0,-1); _stencil_val(u.x,0,0,-1);
#line 229
_stencil_val_a(taux.x,0,0,0);    
        
      }}end__stencil_is_face_x()end_foreach_face_stencil();{
#line 228
foreach_face_generic()is_face_x(){
 val(taux.x,0,0,0) = val(mu.x,0,0,0)*(val(u.z,0,0,0) - val(u.z,-1,0,0) +
      (val(u.x,-1,0,1) + val(u.x,0,0,1))/4. -
      (val(u.x,-1,0,-1) + val(u.x,0,0,-1))/4.)/Delta;}end_is_face_x()end_foreach_face_generic();}}else {_coord _const_mu={_constant[mu.z.i-_NVARMAX],_constant[mu.x.i-_NVARMAX],_constant[mu.y.i-_NVARMAX]};NOT_UNUSED(_const_mu);

      
#line 228
foreach_face_stencil(1,{(NonLocal[]){{"u","vector",(void *)&u,NULL,0},{"_const_mu","_coord",(void *)&_const_mu,NULL,0},{"taux","vector",(void *)&taux,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n// #line 228 \"/home/lisergey/basilisk/src/viscosity.h\"\nis_face_x(){\n val_out_(taux.x,0,0,0) = _const_mu.x*(val(u.z,0,0,0) - val(u.z,-1,0,0) +\n      (val(u.x,-1,0,1) + val(u.x,0,0,1))/4. -\n      (val(u.x,-1,0,-1) + val(u.x,0,0,-1))/4.)/Delta;}end_is_face_x()"})_stencil_is_face_x(){
 {;_stencil_val(u.z,0,0,0); _stencil_val(u.z,-1,0,0);
_stencil_val(u.x,-1,0,1); _stencil_val(u.x,0,0,1);
_stencil_val(u.x,-1,0,-1); _stencil_val(u.x,0,0,-1);
#line 229
_stencil_val_a(taux.x,0,0,0);    
        
      }}end__stencil_is_face_x()end_foreach_face_stencil();

      {
#line 228
foreach_face_generic()is_face_x(){
 val(taux.x,0,0,0) = _const_mu.x*(val(u.z,0,0,0) - val(u.z,-1,0,0) +
      (val(u.x,-1,0,1) + val(u.x,0,0,1))/4. -
      (val(u.x,-1,0,-1) + val(u.x,0,0,-1))/4.)/Delta;}end_is_face_x()end_foreach_face_generic();}}


      if(!is_constant(mu.z)){


      
#line 234
foreach_face_stencil(1,{(NonLocal[]){{"u","vector",(void *)&u,NULL,0},{"mu","vector",(void *)&mu,NULL,0},{"taux","vector",(void *)&taux,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n// #line 234 \"/home/lisergey/basilisk/src/viscosity.h\"\nis_face_y(){\n val_out_(taux.y,0,0,0) = val(mu.y,0,0,0)*(val(u.z,0,0,0) - val(u.z,0,-1,0) +\n      (val(u.y,0,-1,1) + val(u.y,0,0,1))/4. -\n      (val(u.y,0,-1,-1) + val(u.y,0,0,-1))/4.)/Delta;}end_is_face_y()"})_stencil_is_face_y(){
 { _stencil_val(mu.y,0,0,0);_stencil_val(u.z,0,0,0); _stencil_val(u.z,0,-1,0);
_stencil_val(u.y,0,-1,1); _stencil_val(u.y,0,0,1);
_stencil_val(u.y,0,-1,-1); _stencil_val(u.y,0,0,-1);
#line 235
_stencil_val_a(taux.y,0,0,0);    
        
      }}end__stencil_is_face_y()end_foreach_face_stencil();{
#line 234
foreach_face_generic()is_face_y(){
 val(taux.y,0,0,0) = val(mu.y,0,0,0)*(val(u.z,0,0,0) - val(u.z,0,-1,0) +
      (val(u.y,0,-1,1) + val(u.y,0,0,1))/4. -
      (val(u.y,0,-1,-1) + val(u.y,0,0,-1))/4.)/Delta;}end_is_face_y()end_foreach_face_generic();}}else {_coord _const_mu={_constant[mu.z.i-_NVARMAX],_constant[mu.x.i-_NVARMAX],_constant[mu.y.i-_NVARMAX]};NOT_UNUSED(_const_mu);


      
#line 234
foreach_face_stencil(1,{(NonLocal[]){{"u","vector",(void *)&u,NULL,0},{"_const_mu","_coord",(void *)&_const_mu,NULL,0},{"taux","vector",(void *)&taux,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n// #line 234 \"/home/lisergey/basilisk/src/viscosity.h\"\nis_face_y(){\n val_out_(taux.y,0,0,0) = _const_mu.y*(val(u.z,0,0,0) - val(u.z,0,-1,0) +\n      (val(u.y,0,-1,1) + val(u.y,0,0,1))/4. -\n      (val(u.y,0,-1,-1) + val(u.y,0,0,-1))/4.)/Delta;}end_is_face_y()"})_stencil_is_face_y(){
 {;_stencil_val(u.z,0,0,0); _stencil_val(u.z,0,-1,0);
_stencil_val(u.y,0,-1,1); _stencil_val(u.y,0,0,1);
_stencil_val(u.y,0,-1,-1); _stencil_val(u.y,0,0,-1);
#line 235
_stencil_val_a(taux.y,0,0,0);    
        
      }}end__stencil_is_face_y()end_foreach_face_stencil();


      {
#line 234
foreach_face_generic()is_face_y(){
 val(taux.y,0,0,0) = _const_mu.y*(val(u.z,0,0,0) - val(u.z,0,-1,0) +
      (val(u.y,0,-1,1) + val(u.y,0,0,1))/4. -
      (val(u.y,0,-1,-1) + val(u.y,0,0,-1))/4.)/Delta;}end_is_face_y()end_foreach_face_generic();}}

    if(!is_constant(rho)){

    
#line 239
foreach_stencil (1,{(NonLocal[]){{"maxres","double",(void *)&maxres,NULL,0,'M'},{"rho","scalar",(void *)&rho,NULL,0},{"dt","double",(void *)&dt,NULL,0},{"u","vector",(void *)&u,NULL,0},{"r","vector",(void *)&r,NULL,0},{"res","vector",(void *)&res,NULL,0},{"taux","vector",(void *)&taux,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},""," \n// #line 239 \"/home/lisergey/basilisk/src/viscosity.h\"\n{\n      real d = 0.;\n      \n d += val(taux.z,0,0,1) - val(taux.z,0,0,0);\n \n// #line 242\nd += val(taux.x,1,0,0) - val(taux.x,0,0,0);\n \n// #line 242\nd += val(taux.y,0,1,0) - val(taux.y,0,0,0);\n      val_out_(res.z,0,0,0) = val(r.z,0,0,0) - ((coord){1.,1.,1.}).z*val(u.z,0,0,0) + dt/val(rho,0,0,0)*d/Delta;\n      if (fabs (val(res.z,0,0,0)) > maxres)\n maxres = fabs (val(res.z,0,0,0));\n    }"}) {   
      
      
 { _stencil_val(taux.z,0,0,1); _stencil_val(taux.z,0,0,0);  }
 
#line 242
{ _stencil_val(taux.x,1,0,0); _stencil_val(taux.x,0,0,0);  }
 
#line 242
{ _stencil_val(taux.y,0,1,0); _stencil_val(taux.y,0,0,0);  } 
_stencil_val(r.z,0,0,0);_stencil_val(u.z,0,0,0);_stencil_val(rho,0,0,0);
      
#line 243
_stencil_val_a(res.z,0,0,0);
_stencil_val(res.z,0,0,0);
 {_stencil_val(res.z,0,0,0);   }     
          
    
#line 246
}end_foreach_stencil();
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel reduction(max:maxres)){
#line 239
foreach () {
      double d = 0.;
      
 d += val(taux.z,0,0,1) - val(taux.z,0,0,0);
 
#line 242
d += val(taux.x,1,0,0) - val(taux.x,0,0,0);
 
#line 242
d += val(taux.y,0,1,0) - val(taux.y,0,0,0);
      val(res.z,0,0,0) = val(r.z,0,0,0) - ((coord){1.,1.,1.}).z*val(u.z,0,0,0) + dt/val(rho,0,0,0)*d/Delta;
      if (fabs (val(res.z,0,0,0)) > maxres)
 maxres = fabs (val(res.z,0,0,0));
    }end_foreach();mpi_all_reduce_array(&maxres,double,MPI_MAX,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 246
}else {double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);

    
#line 239
foreach_stencil (1,{(NonLocal[]){{"maxres","double",(void *)&maxres,NULL,0,'M'},{"_const_rho","double",(void *)&_const_rho,NULL,0},{"dt","double",(void *)&dt,NULL,0},{"u","vector",(void *)&u,NULL,0},{"r","vector",(void *)&r,NULL,0},{"res","vector",(void *)&res,NULL,0},{"taux","vector",(void *)&taux,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},""," \n// #line 239 \"/home/lisergey/basilisk/src/viscosity.h\"\n{\n      real d = 0.;\n      \n d += val(taux.z,0,0,1) - val(taux.z,0,0,0);\n \n// #line 242\nd += val(taux.x,1,0,0) - val(taux.x,0,0,0);\n \n// #line 242\nd += val(taux.y,0,1,0) - val(taux.y,0,0,0);\n      val_out_(res.z,0,0,0) = val(r.z,0,0,0) - ((coord){1.,1.,1.}).z*val(u.z,0,0,0) + dt/_const_rho*d/Delta;\n      if (fabs (val(res.z,0,0,0)) > maxres)\n maxres = fabs (val(res.z,0,0,0));\n    }"}) {   
      
      
 { _stencil_val(taux.z,0,0,1); _stencil_val(taux.z,0,0,0);  }
 
#line 242
{ _stencil_val(taux.x,1,0,0); _stencil_val(taux.x,0,0,0);  }
 
#line 242
{ _stencil_val(taux.y,0,1,0); _stencil_val(taux.y,0,0,0);  } 
_stencil_val(r.z,0,0,0);_stencil_val(u.z,0,0,0);;
      
#line 243
_stencil_val_a(res.z,0,0,0);
_stencil_val(res.z,0,0,0);
 {_stencil_val(res.z,0,0,0);   }     
          
    
#line 246
}end_foreach_stencil();

    
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel reduction(max:maxres)){
#line 239
foreach () {
      double d = 0.;
      
 d += val(taux.z,0,0,1) - val(taux.z,0,0,0);
 
#line 242
d += val(taux.x,1,0,0) - val(taux.x,0,0,0);
 
#line 242
d += val(taux.y,0,1,0) - val(taux.y,0,0,0);
      val(res.z,0,0,0) = val(r.z,0,0,0) - ((coord){1.,1.,1.}).z*val(u.z,0,0,0) + dt/_const_rho*d/Delta;
      if (fabs (val(res.z,0,0,0)) > maxres)
 maxres = fabs (val(res.z,0,0,0));
    }end_foreach();mpi_all_reduce_array(&maxres,double,MPI_MAX,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 246
}delete((scalar*)((vector[]){taux,{{-1},{-1},{-1}}}));
  }
#line 276 "/home/lisergey/basilisk/src/viscosity.h"
  return maxres;
}
#line 288 "/home/lisergey/basilisk/src/viscosity.h"
     
mgstats viscosity (vector u, vector mu, scalar rho, double dt,
     int nrelax, scalar * res)
{tracing("viscosity","/home/lisergey/basilisk/src/viscosity.h",289);





  vector  r=new_vector("r");
  foreach_stencil(1,{(NonLocal[]){{"u","vector",(void *)&u,NULL,0},{"r","vector",(void *)&r,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n    \n// #line 299 \"/home/lisergey/basilisk/src/viscosity.h\"\n{\n      val_out_(r.x,0,0,0) = val(u.x,0,0,0);\n      \n// #line 300\nval_out_(r.y,0,0,0) = val(u.y,0,0,0);\n      \n// #line 300\nval_out_(r.z,0,0,0) = val(u.z,0,0,0);}"})
    {
      { _stencil_val(u.x,0,0,0);_stencil_val_a(r.x,0,0,0); }
      
#line 300
{ _stencil_val(u.y,0,0,0);_stencil_val_a(r.y,0,0,0); }
      
#line 300
{ _stencil_val(u.z,0,0,0);_stencil_val_a(r.z,0,0,0); }}end_foreach_stencil();
  {
#line 298
foreach()
    {
      val(r.x,0,0,0) = val(u.x,0,0,0);
      
#line 300
val(r.y,0,0,0) = val(u.y,0,0,0);
      
#line 300
val(r.z,0,0,0) = val(u.z,0,0,0);}end_foreach();}




  restriction (((scalar[]){mu.x,mu.y,mu.z,rho,{-1}}));
  struct Viscosity p = { mu, rho, dt };
  { mgstats _ret= mg_solve ((scalar *)((vector[]){u,{{-1},{-1},{-1}}}), (scalar *)((vector[]){r,{{-1},{-1},{-1}}})
,
     
#line 308
residual_viscosity, relax_viscosity, &p, nrelax, res
#line 132 "/home/lisergey/basilisk/src/poisson.h"
, 
0, 
TOLERANCE
#line 308 "/home/lisergey/basilisk/src/viscosity.h"
);delete((scalar*)((vector[]){r,{{-1},{-1},{-1}}}));{end_tracing("viscosity","/home/lisergey/basilisk/src/viscosity.h",308);return _ret;}}delete((scalar*)((vector[]){r,{{-1},{-1},{-1}}}));
end_tracing("viscosity","/home/lisergey/basilisk/src/viscosity.h",309);}
#line 318 "/home/lisergey/basilisk/src/viscosity.h"
     
mgstats viscosity_explicit (vector u, vector mu, scalar rho, double dt)
{tracing("viscosity_explicit","/home/lisergey/basilisk/src/viscosity.h",319);
  vector  r=new_vector("r");
  mgstats mg = {0};
  struct Viscosity p = { mu, rho, dt };
  mg.resb = residual_viscosity ((scalar *)((vector[]){u,{{-1},{-1},{-1}}}), (scalar *)((vector[]){u,{{-1},{-1},{-1}}}), (scalar *)((vector[]){r,{{-1},{-1},{-1}}}), &p);
  foreach_stencil(1,{(NonLocal[]){{"r","vector",(void *)&r,NULL,0},{"u","vector",(void *)&u,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n    \n// #line 326 \"/home/lisergey/basilisk/src/viscosity.h\"\n{\n      val_out_(u.x,0,0,0) += val(r.x,0,0,0);\n      \n// #line 327\nval_out_(u.y,0,0,0) += val(r.y,0,0,0);\n      \n// #line 327\nval_out_(u.z,0,0,0) += val(r.z,0,0,0);}"})
    {
      { _stencil_val(r.x,0,0,0);_stencil_val_r(u.x,0,0,0); }
      
#line 327
{ _stencil_val(r.y,0,0,0);_stencil_val_r(u.y,0,0,0); }
      
#line 327
{ _stencil_val(r.z,0,0,0);_stencil_val_r(u.z,0,0,0); }}end_foreach_stencil();
  {
#line 325
foreach()
    {
      val(u.x,0,0,0) += val(r.x,0,0,0);
      
#line 327
val(u.y,0,0,0) += val(r.y,0,0,0);
      
#line 327
val(u.z,0,0,0) += val(r.z,0,0,0);}end_foreach();}
  {delete((scalar*)((vector[]){r,{{-1},{-1},{-1}}}));{end_tracing("viscosity_explicit","/home/lisergey/basilisk/src/viscosity.h",328);return mg;}}delete((scalar*)((vector[]){r,{{-1},{-1},{-1}}}));
end_tracing("viscosity_explicit","/home/lisergey/basilisk/src/viscosity.h",329);}
#line 34 "/home/lisergey/basilisk/src/navier-stokes/centered.h"
#line 44 "/home/lisergey/basilisk/src/navier-stokes/centered.h"
scalar  p={0};
vector  u={{1},{2},{3}},  g={{4},{5},{6}};
scalar  pf={7};
vector  uf={{8},{9},{10}};
#line 70 "/home/lisergey/basilisk/src/navier-stokes/centered.h"
        vector mu = {{_NVARMAX+0},{_NVARMAX+1},{_NVARMAX+2}}, a = {{_NVARMAX+0},{_NVARMAX+1},{_NVARMAX+2}}, alpha = {{_NVARMAX+3},{_NVARMAX+4},{_NVARMAX+5}};
        scalar rho = {_NVARMAX+6};
mgstats mgp = {0}, mgpf = {0}, mgu = {0};
bool stokes = false;
#line 91
static double _boundary0(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;if(!is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)){{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann((val(a.x,1,0,0)*val(fm.x,1,0,0)/val(alpha.x,1,0,0)), point, neighbor, _s, data));}}}}else if(is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX],_constant[a.z.i-_NVARMAX]};NOT_UNUSED(_const_a);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann((_const_a.x*val(fm.x,1,0,0)/val(alpha.x,1,0,0)), point, neighbor, _s, data));}}}}else if(!is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann((val(a.x,1,0,0)*_const_fm.x/val(alpha.x,1,0,0)), point, neighbor, _s, data));}}}}else if(is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX],_constant[a.z.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann((_const_a.x*_const_fm.x/val(alpha.x,1,0,0)), point, neighbor, _s, data));}}}}else if(!is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX],_constant[alpha.z.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann((val(a.x,1,0,0)*val(fm.x,1,0,0)/_const_alpha.x), point, neighbor, _s, data));}}}}else if(is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX],_constant[a.z.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX],_constant[alpha.z.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann((_const_a.x*val(fm.x,1,0,0)/_const_alpha.x), point, neighbor, _s, data));}}}}else if(!is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX],_constant[alpha.z.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann((val(a.x,1,0,0)*_const_fm.x/_const_alpha.x), point, neighbor, _s, data));}}}}else {_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX],_constant[a.z.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX],_constant[alpha.z.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann((_const_a.x*_const_fm.x/_const_alpha.x), point, neighbor, _s, data));}}}}}static double _boundary0_homogeneous(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;if(!is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)){{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann_homogeneous((val(a.x,1,0,0)*val(fm.x,1,0,0)/val(alpha.x,1,0,0)), point, neighbor, _s, data));}}}}else if(is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX],_constant[a.z.i-_NVARMAX]};NOT_UNUSED(_const_a);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann_homogeneous((_const_a.x*val(fm.x,1,0,0)/val(alpha.x,1,0,0)), point, neighbor, _s, data));}}}}else if(!is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann_homogeneous((val(a.x,1,0,0)*_const_fm.x/val(alpha.x,1,0,0)), point, neighbor, _s, data));}}}}else if(is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX],_constant[a.z.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann_homogeneous((_const_a.x*_const_fm.x/val(alpha.x,1,0,0)), point, neighbor, _s, data));}}}}else if(!is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX],_constant[alpha.z.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann_homogeneous((val(a.x,1,0,0)*val(fm.x,1,0,0)/_const_alpha.x), point, neighbor, _s, data));}}}}else if(is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX],_constant[a.z.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX],_constant[alpha.z.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann_homogeneous((_const_a.x*val(fm.x,1,0,0)/_const_alpha.x), point, neighbor, _s, data));}}}}else if(!is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX],_constant[alpha.z.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann_homogeneous((val(a.x,1,0,0)*_const_fm.x/_const_alpha.x), point, neighbor, _s, data));}}}}else {_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX],_constant[a.z.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX],_constant[alpha.z.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann_homogeneous((_const_a.x*_const_fm.x/_const_alpha.x), point, neighbor, _s, data));}}}}}
static double _boundary1(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;if(!is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)){{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann(- (val(a.x,0,0,0)*val(fm.x,0,0,0)/val(alpha.x,0,0,0)), point, neighbor, _s, data));}}}}else if(is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX],_constant[a.z.i-_NVARMAX]};NOT_UNUSED(_const_a);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann(- (_const_a.x*val(fm.x,0,0,0)/val(alpha.x,0,0,0)), point, neighbor, _s, data));}}}}else if(!is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann(- (val(a.x,0,0,0)*_const_fm.x/val(alpha.x,0,0,0)), point, neighbor, _s, data));}}}}else if(is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX],_constant[a.z.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann(- (_const_a.x*_const_fm.x/val(alpha.x,0,0,0)), point, neighbor, _s, data));}}}}else if(!is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX],_constant[alpha.z.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann(- (val(a.x,0,0,0)*val(fm.x,0,0,0)/_const_alpha.x), point, neighbor, _s, data));}}}}else if(is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX],_constant[a.z.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX],_constant[alpha.z.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann(- (_const_a.x*val(fm.x,0,0,0)/_const_alpha.x), point, neighbor, _s, data));}}}}else if(!is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX],_constant[alpha.z.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann(- (val(a.x,0,0,0)*_const_fm.x/_const_alpha.x), point, neighbor, _s, data));}}}}else {_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX],_constant[a.z.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX],_constant[alpha.z.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann(- (_const_a.x*_const_fm.x/_const_alpha.x), point, neighbor, _s, data));}}}}}static double _boundary1_homogeneous(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;if(!is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)){{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann_homogeneous(- (val(a.x,0,0,0)*val(fm.x,0,0,0)/val(alpha.x,0,0,0)), point, neighbor, _s, data));}}}}else if(is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX],_constant[a.z.i-_NVARMAX]};NOT_UNUSED(_const_a);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann_homogeneous(- (_const_a.x*val(fm.x,0,0,0)/val(alpha.x,0,0,0)), point, neighbor, _s, data));}}}}else if(!is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann_homogeneous(- (val(a.x,0,0,0)*_const_fm.x/val(alpha.x,0,0,0)), point, neighbor, _s, data));}}}}else if(is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX],_constant[a.z.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann_homogeneous(- (_const_a.x*_const_fm.x/val(alpha.x,0,0,0)), point, neighbor, _s, data));}}}}else if(!is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX],_constant[alpha.z.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann_homogeneous(- (val(a.x,0,0,0)*val(fm.x,0,0,0)/_const_alpha.x), point, neighbor, _s, data));}}}}else if(is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX],_constant[a.z.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX],_constant[alpha.z.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann_homogeneous(- (_const_a.x*val(fm.x,0,0,0)/_const_alpha.x), point, neighbor, _s, data));}}}}else if(!is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX],_constant[alpha.z.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann_homogeneous(- (val(a.x,0,0,0)*_const_fm.x/_const_alpha.x), point, neighbor, _s, data));}}}}else {_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX],_constant[a.z.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX],_constant[alpha.z.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann_homogeneous(- (_const_a.x*_const_fm.x/_const_alpha.x), point, neighbor, _s, data));}}}}}








static double _boundary2(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;if(!is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)){{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann((val(a.y,0,1,0)*val(fm.y,0,1,0)/val(alpha.y,0,1,0)), point, neighbor, _s, data));}}}}else if(is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX],_constant[a.z.i-_NVARMAX]};NOT_UNUSED(_const_a);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann((_const_a.y*val(fm.y,0,1,0)/val(alpha.y,0,1,0)), point, neighbor, _s, data));}}}}else if(!is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann((val(a.y,0,1,0)*_const_fm.y/val(alpha.y,0,1,0)), point, neighbor, _s, data));}}}}else if(is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX],_constant[a.z.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann((_const_a.y*_const_fm.y/val(alpha.y,0,1,0)), point, neighbor, _s, data));}}}}else if(!is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX],_constant[alpha.z.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann((val(a.y,0,1,0)*val(fm.y,0,1,0)/_const_alpha.y), point, neighbor, _s, data));}}}}else if(is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX],_constant[a.z.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX],_constant[alpha.z.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann((_const_a.y*val(fm.y,0,1,0)/_const_alpha.y), point, neighbor, _s, data));}}}}else if(!is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX],_constant[alpha.z.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann((val(a.y,0,1,0)*_const_fm.y/_const_alpha.y), point, neighbor, _s, data));}}}}else {_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX],_constant[a.z.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX],_constant[alpha.z.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann((_const_a.y*_const_fm.y/_const_alpha.y), point, neighbor, _s, data));}}}}}static double _boundary2_homogeneous(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;if(!is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)){{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann_homogeneous((val(a.y,0,1,0)*val(fm.y,0,1,0)/val(alpha.y,0,1,0)), point, neighbor, _s, data));}}}}else if(is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX],_constant[a.z.i-_NVARMAX]};NOT_UNUSED(_const_a);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann_homogeneous((_const_a.y*val(fm.y,0,1,0)/val(alpha.y,0,1,0)), point, neighbor, _s, data));}}}}else if(!is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann_homogeneous((val(a.y,0,1,0)*_const_fm.y/val(alpha.y,0,1,0)), point, neighbor, _s, data));}}}}else if(is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX],_constant[a.z.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann_homogeneous((_const_a.y*_const_fm.y/val(alpha.y,0,1,0)), point, neighbor, _s, data));}}}}else if(!is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX],_constant[alpha.z.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann_homogeneous((val(a.y,0,1,0)*val(fm.y,0,1,0)/_const_alpha.y), point, neighbor, _s, data));}}}}else if(is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX],_constant[a.z.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX],_constant[alpha.z.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann_homogeneous((_const_a.y*val(fm.y,0,1,0)/_const_alpha.y), point, neighbor, _s, data));}}}}else if(!is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX],_constant[alpha.z.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann_homogeneous((val(a.y,0,1,0)*_const_fm.y/_const_alpha.y), point, neighbor, _s, data));}}}}else {_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX],_constant[a.z.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX],_constant[alpha.z.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann_homogeneous((_const_a.y*_const_fm.y/_const_alpha.y), point, neighbor, _s, data));}}}}}
static double _boundary3(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;if(!is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)){{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann(- (val(a.y,0,0,0)*val(fm.y,0,0,0)/val(alpha.y,0,0,0)), point, neighbor, _s, data));}}}}else if(is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX],_constant[a.z.i-_NVARMAX]};NOT_UNUSED(_const_a);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann(- (_const_a.y*val(fm.y,0,0,0)/val(alpha.y,0,0,0)), point, neighbor, _s, data));}}}}else if(!is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann(- (val(a.y,0,0,0)*_const_fm.y/val(alpha.y,0,0,0)), point, neighbor, _s, data));}}}}else if(is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX],_constant[a.z.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann(- (_const_a.y*_const_fm.y/val(alpha.y,0,0,0)), point, neighbor, _s, data));}}}}else if(!is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX],_constant[alpha.z.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann(- (val(a.y,0,0,0)*val(fm.y,0,0,0)/_const_alpha.y), point, neighbor, _s, data));}}}}else if(is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX],_constant[a.z.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX],_constant[alpha.z.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann(- (_const_a.y*val(fm.y,0,0,0)/_const_alpha.y), point, neighbor, _s, data));}}}}else if(!is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX],_constant[alpha.z.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann(- (val(a.y,0,0,0)*_const_fm.y/_const_alpha.y), point, neighbor, _s, data));}}}}else {_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX],_constant[a.z.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX],_constant[alpha.z.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann(- (_const_a.y*_const_fm.y/_const_alpha.y), point, neighbor, _s, data));}}}}}static double _boundary3_homogeneous(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;if(!is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)){{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann_homogeneous(- (val(a.y,0,0,0)*val(fm.y,0,0,0)/val(alpha.y,0,0,0)), point, neighbor, _s, data));}}}}else if(is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX],_constant[a.z.i-_NVARMAX]};NOT_UNUSED(_const_a);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann_homogeneous(- (_const_a.y*val(fm.y,0,0,0)/val(alpha.y,0,0,0)), point, neighbor, _s, data));}}}}else if(!is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann_homogeneous(- (val(a.y,0,0,0)*_const_fm.y/val(alpha.y,0,0,0)), point, neighbor, _s, data));}}}}else if(is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX],_constant[a.z.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann_homogeneous(- (_const_a.y*_const_fm.y/val(alpha.y,0,0,0)), point, neighbor, _s, data));}}}}else if(!is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX],_constant[alpha.z.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann_homogeneous(- (val(a.y,0,0,0)*val(fm.y,0,0,0)/_const_alpha.y), point, neighbor, _s, data));}}}}else if(is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX],_constant[a.z.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX],_constant[alpha.z.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann_homogeneous(- (_const_a.y*val(fm.y,0,0,0)/_const_alpha.y), point, neighbor, _s, data));}}}}else if(!is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX],_constant[alpha.z.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann_homogeneous(- (val(a.y,0,0,0)*_const_fm.y/_const_alpha.y), point, neighbor, _s, data));}}}}else {_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX],_constant[a.z.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX],_constant[alpha.z.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann_homogeneous(- (_const_a.y*_const_fm.y/_const_alpha.y), point, neighbor, _s, data));}}}}}


static double _boundary4(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;if(!is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)){{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann((val(a.z,0,0,1)*val(fm.z,0,0,1)/val(alpha.z,0,0,1)), point, neighbor, _s, data));}}}}else if(is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX],_constant[a.z.i-_NVARMAX]};NOT_UNUSED(_const_a);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann((_const_a.z*val(fm.z,0,0,1)/val(alpha.z,0,0,1)), point, neighbor, _s, data));}}}}else if(!is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann((val(a.z,0,0,1)*_const_fm.z/val(alpha.z,0,0,1)), point, neighbor, _s, data));}}}}else if(is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX],_constant[a.z.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann((_const_a.z*_const_fm.z/val(alpha.z,0,0,1)), point, neighbor, _s, data));}}}}else if(!is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX],_constant[alpha.z.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann((val(a.z,0,0,1)*val(fm.z,0,0,1)/_const_alpha.z), point, neighbor, _s, data));}}}}else if(is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX],_constant[a.z.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX],_constant[alpha.z.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann((_const_a.z*val(fm.z,0,0,1)/_const_alpha.z), point, neighbor, _s, data));}}}}else if(!is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX],_constant[alpha.z.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann((val(a.z,0,0,1)*_const_fm.z/_const_alpha.z), point, neighbor, _s, data));}}}}else {_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX],_constant[a.z.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX],_constant[alpha.z.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann((_const_a.z*_const_fm.z/_const_alpha.z), point, neighbor, _s, data));}}}}}static double _boundary4_homogeneous(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;if(!is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)){{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann_homogeneous((val(a.z,0,0,1)*val(fm.z,0,0,1)/val(alpha.z,0,0,1)), point, neighbor, _s, data));}}}}else if(is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX],_constant[a.z.i-_NVARMAX]};NOT_UNUSED(_const_a);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann_homogeneous((_const_a.z*val(fm.z,0,0,1)/val(alpha.z,0,0,1)), point, neighbor, _s, data));}}}}else if(!is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann_homogeneous((val(a.z,0,0,1)*_const_fm.z/val(alpha.z,0,0,1)), point, neighbor, _s, data));}}}}else if(is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX],_constant[a.z.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann_homogeneous((_const_a.z*_const_fm.z/val(alpha.z,0,0,1)), point, neighbor, _s, data));}}}}else if(!is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX],_constant[alpha.z.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann_homogeneous((val(a.z,0,0,1)*val(fm.z,0,0,1)/_const_alpha.z), point, neighbor, _s, data));}}}}else if(is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX],_constant[a.z.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX],_constant[alpha.z.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann_homogeneous((_const_a.z*val(fm.z,0,0,1)/_const_alpha.z), point, neighbor, _s, data));}}}}else if(!is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX],_constant[alpha.z.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann_homogeneous((val(a.z,0,0,1)*_const_fm.z/_const_alpha.z), point, neighbor, _s, data));}}}}else {_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX],_constant[a.z.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX],_constant[alpha.z.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann_homogeneous((_const_a.z*_const_fm.z/_const_alpha.z), point, neighbor, _s, data));}}}}}
static double _boundary5(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;if(!is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)){{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann(- (val(a.z,0,0,0)*val(fm.z,0,0,0)/val(alpha.z,0,0,0)), point, neighbor, _s, data));}}}}else if(is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX],_constant[a.z.i-_NVARMAX]};NOT_UNUSED(_const_a);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann(- (_const_a.z*val(fm.z,0,0,0)/val(alpha.z,0,0,0)), point, neighbor, _s, data));}}}}else if(!is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann(- (val(a.z,0,0,0)*_const_fm.z/val(alpha.z,0,0,0)), point, neighbor, _s, data));}}}}else if(is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX],_constant[a.z.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann(- (_const_a.z*_const_fm.z/val(alpha.z,0,0,0)), point, neighbor, _s, data));}}}}else if(!is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX],_constant[alpha.z.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann(- (val(a.z,0,0,0)*val(fm.z,0,0,0)/_const_alpha.z), point, neighbor, _s, data));}}}}else if(is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX],_constant[a.z.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX],_constant[alpha.z.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann(- (_const_a.z*val(fm.z,0,0,0)/_const_alpha.z), point, neighbor, _s, data));}}}}else if(!is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX],_constant[alpha.z.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann(- (val(a.z,0,0,0)*_const_fm.z/_const_alpha.z), point, neighbor, _s, data));}}}}else {_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX],_constant[a.z.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX],_constant[alpha.z.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann(- (_const_a.z*_const_fm.z/_const_alpha.z), point, neighbor, _s, data));}}}}}static double _boundary5_homogeneous(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;if(!is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)){{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann_homogeneous(- (val(a.z,0,0,0)*val(fm.z,0,0,0)/val(alpha.z,0,0,0)), point, neighbor, _s, data));}}}}else if(is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX],_constant[a.z.i-_NVARMAX]};NOT_UNUSED(_const_a);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann_homogeneous(- (_const_a.z*val(fm.z,0,0,0)/val(alpha.z,0,0,0)), point, neighbor, _s, data));}}}}else if(!is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann_homogeneous(- (val(a.z,0,0,0)*_const_fm.z/val(alpha.z,0,0,0)), point, neighbor, _s, data));}}}}else if(is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX],_constant[a.z.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann_homogeneous(- (_const_a.z*_const_fm.z/val(alpha.z,0,0,0)), point, neighbor, _s, data));}}}}else if(!is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX],_constant[alpha.z.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann_homogeneous(- (val(a.z,0,0,0)*val(fm.z,0,0,0)/_const_alpha.z), point, neighbor, _s, data));}}}}else if(is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX],_constant[a.z.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX],_constant[alpha.z.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann_homogeneous(- (_const_a.z*val(fm.z,0,0,0)/_const_alpha.z), point, neighbor, _s, data));}}}}else if(!is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX],_constant[alpha.z.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann_homogeneous(- (val(a.z,0,0,0)*_const_fm.z/_const_alpha.z), point, neighbor, _s, data));}}}}else {_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX],_constant[a.z.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX],_constant[alpha.z.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann_homogeneous(- (_const_a.z*_const_fm.z/_const_alpha.z), point, neighbor, _s, data));}}}}}
#line 126
static int defaults_0_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i = 0)!=0;*ip=i;*tp=t;return ret;}
#line 126 "/home/lisergey/basilisk/src/navier-stokes/centered.h"
      static int defaults_0(const int i,const double t,Event *_ev){tracing("defaults_0","/home/lisergey/basilisk/src/navier-stokes/centered.h",126);
{

  CFL = 0.8;




  _attribute[p.i].nodump = _attribute[pf.i].nodump = true;





  if (alpha.x.i == unityf.x.i) {
    alpha = fm;
    rho = cm;
  }
  else if (!is_constant(alpha.x)) {
    vector alphav = alpha;
    if(!is_constant(fm.x)){
    
#line 146
foreach_face_stencil(1,{(NonLocal[]){{"fm","vector",(void *)&fm,NULL,0},{"alphav","vector",(void *)&alphav,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n// #line 146 \"/home/lisergey/basilisk/src/navier-stokes/centered.h\"\n{is_face_x(){\n      val_out_(alphav.x,0,0,0) = val(fm.x,0,0,0);}end_is_face_x()\n// #line 146\nis_face_y(){\n      val_out_(alphav.y,0,0,0) = val(fm.y,0,0,0);}end_is_face_y()\n// #line 146\nis_face_z(){\n      val_out_(alphav.z,0,0,0) = val(fm.z,0,0,0);}end_is_face_z()}"}){_stencil_is_face_x(){
      { _stencil_val(fm.x,0,0,0);_stencil_val_a(alphav.x,0,0,0); }}end__stencil_is_face_x()
#line 146
_stencil_is_face_y(){
      { _stencil_val(fm.y,0,0,0);_stencil_val_a(alphav.y,0,0,0); }}end__stencil_is_face_y()
#line 146
_stencil_is_face_z(){
      { _stencil_val(fm.z,0,0,0);_stencil_val_a(alphav.z,0,0,0); }}end__stencil_is_face_z()}end_foreach_face_stencil();{
#line 146
foreach_face_generic(){is_face_x(){
      val(alphav.x,0,0,0) = val(fm.x,0,0,0);}end_is_face_x()
#line 146
is_face_y(){
      val(alphav.y,0,0,0) = val(fm.y,0,0,0);}end_is_face_y()
#line 146
is_face_z(){
      val(alphav.z,0,0,0) = val(fm.z,0,0,0);}end_is_face_z()}end_foreach_face_generic();}}else {_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);
    
#line 146
foreach_face_stencil(1,{(NonLocal[]){{"_const_fm","_coord",(void *)&_const_fm,NULL,0},{"alphav","vector",(void *)&alphav,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n// #line 146 \"/home/lisergey/basilisk/src/navier-stokes/centered.h\"\n{is_face_x(){\n      val_out_(alphav.x,0,0,0) = _const_fm.x;}end_is_face_x()\n// #line 146\nis_face_y(){\n      val_out_(alphav.y,0,0,0) = _const_fm.y;}end_is_face_y()\n// #line 146\nis_face_z(){\n      val_out_(alphav.z,0,0,0) = _const_fm.z;}end_is_face_z()}"}){_stencil_is_face_x(){
      {;_stencil_val_a(alphav.x,0,0,0); }}end__stencil_is_face_x()
#line 146
_stencil_is_face_y(){
      {;_stencil_val_a(alphav.y,0,0,0); }}end__stencil_is_face_y()
#line 146
_stencil_is_face_z(){
      {;_stencil_val_a(alphav.z,0,0,0); }}end__stencil_is_face_z()}end_foreach_face_stencil();
    {
#line 146
foreach_face_generic(){is_face_x(){
      val(alphav.x,0,0,0) = _const_fm.x;}end_is_face_x()
#line 146
is_face_y(){
      val(alphav.y,0,0,0) = _const_fm.y;}end_is_face_y()
#line 146
is_face_z(){
      val(alphav.z,0,0,0) = _const_fm.z;}end_is_face_z()}end_foreach_face_generic();}}
  }






  _attribute[uf.x.i].refine = refine_face_solenoidal;
#line 178 "/home/lisergey/basilisk/src/navier-stokes/centered.h"
  foreach_stencil(1,{(NonLocal[]){{"t","double",(void *)&t,NULL,0},{"u","vector",(void *)&u,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n    \n// #line 179 \"/home/lisergey/basilisk/src/navier-stokes/centered.h\"\n{\n      dimensional (val(u.x,0,0,0) == Delta/t);\n      \n// #line 180\ndimensional (val(u.y,0,0,0) == Delta/t);\n      \n// #line 180\ndimensional (val(u.z,0,0,0) == Delta/t);}"})
    {
      {_stencil_val(u.x,0,0,0);   }
      
#line 180
{_stencil_val(u.y,0,0,0);   }
      
#line 180
{_stencil_val(u.z,0,0,0);   }}end_foreach_stencil();
#line 178 "/home/lisergey/basilisk/src/navier-stokes/centered.h"
  {foreach()
    {
      dimensional (val(u.x,0,0,0) == Delta/t);
      
#line 180
dimensional (val(u.y,0,0,0) == Delta/t);
      
#line 180
dimensional (val(u.z,0,0,0) == Delta/t);}end_foreach();}
}{end_tracing("defaults_0","/home/lisergey/basilisk/src/navier-stokes/centered.h",181);return 0;}end_tracing("defaults_0","/home/lisergey/basilisk/src/navier-stokes/centered.h",181);}





static int default_display_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i = 0)!=0;*ip=i;*tp=t;return ret;}






#line 187
      static int default_display(const int i,const double t,Event *_ev){tracing("default_display","/home/lisergey/basilisk/src/navier-stokes/centered.h",187);
  display ("squares (color = 'u.x', spread = -1);"
#line 1425 "/home/lisergey/basilisk/src/common.h"
, false
#line 188 "/home/lisergey/basilisk/src/navier-stokes/centered.h"
);{end_tracing("default_display","/home/lisergey/basilisk/src/navier-stokes/centered.h",188);return 0;}end_tracing("default_display","/home/lisergey/basilisk/src/navier-stokes/centered.h",188);}





double dtmax;

static int init_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i = 0)!=0;*ip=i;*tp=t;return ret;}


#line 196
      static int init(const int i,const double t,Event *_ev){tracing("init","/home/lisergey/basilisk/src/navier-stokes/centered.h",196);
{
  trash (((vector[]){uf,{{-1},{-1},{-1}}}));
  if(!is_constant(fm.x)){
  
#line 199
foreach_face_stencil(1,{(NonLocal[]){{"u","vector",(void *)&u,NULL,0},{"fm","vector",(void *)&fm,NULL,0},{"uf","vector",(void *)&uf,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n// #line 199 \"/home/lisergey/basilisk/src/navier-stokes/centered.h\"\n{is_face_x(){\n    val_out_(uf.x,0,0,0) = val(fm.x,0,0,0)*((val(u.x,0,0,0) + val(u.x,0 -1,0,0))/2.);}end_is_face_x()\n// #line 199\nis_face_y(){\n    val_out_(uf.y,0,0,0) = val(fm.y,0,0,0)*((val(u.y,0,0,0) + val(u.y,0,0 -1,0))/2.);}end_is_face_y()\n// #line 199\nis_face_z(){\n    val_out_(uf.z,0,0,0) = val(fm.z,0,0,0)*((val(u.z,0,0,0) + val(u.z,0,0,0 -1))/2.);}end_is_face_z()}"}){_stencil_is_face_x(){
    { _stencil_val(fm.x,0,0,0);_stencil_val(u.x,0,0,0); _stencil_val(u.x,0 -1,0,0);_stencil_val_a(uf.x,0,0,0);  }}end__stencil_is_face_x()
#line 199
_stencil_is_face_y(){
    { _stencil_val(fm.y,0,0,0);_stencil_val(u.y,0,0,0); _stencil_val(u.y,0,0 -1,0);_stencil_val_a(uf.y,0,0,0);  }}end__stencil_is_face_y()
#line 199
_stencil_is_face_z(){
    { _stencil_val(fm.z,0,0,0);_stencil_val(u.z,0,0,0); _stencil_val(u.z,0,0,0 -1);_stencil_val_a(uf.z,0,0,0);  }}end__stencil_is_face_z()}end_foreach_face_stencil();{
#line 199
foreach_face_generic(){is_face_x(){
    val(uf.x,0,0,0) = val(fm.x,0,0,0)*((val(u.x,0,0,0) + val(u.x,0 -1,0,0))/2.);}end_is_face_x()
#line 199
is_face_y(){
    val(uf.y,0,0,0) = val(fm.y,0,0,0)*((val(u.y,0,0,0) + val(u.y,0,0 -1,0))/2.);}end_is_face_y()
#line 199
is_face_z(){
    val(uf.z,0,0,0) = val(fm.z,0,0,0)*((val(u.z,0,0,0) + val(u.z,0,0,0 -1))/2.);}end_is_face_z()}end_foreach_face_generic();}}else {_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);
  
#line 199
foreach_face_stencil(1,{(NonLocal[]){{"u","vector",(void *)&u,NULL,0},{"_const_fm","_coord",(void *)&_const_fm,NULL,0},{"uf","vector",(void *)&uf,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n// #line 199 \"/home/lisergey/basilisk/src/navier-stokes/centered.h\"\n{is_face_x(){\n    val_out_(uf.x,0,0,0) = _const_fm.x*((val(u.x,0,0,0) + val(u.x,0 -1,0,0))/2.);}end_is_face_x()\n// #line 199\nis_face_y(){\n    val_out_(uf.y,0,0,0) = _const_fm.y*((val(u.y,0,0,0) + val(u.y,0,0 -1,0))/2.);}end_is_face_y()\n// #line 199\nis_face_z(){\n    val_out_(uf.z,0,0,0) = _const_fm.z*((val(u.z,0,0,0) + val(u.z,0,0,0 -1))/2.);}end_is_face_z()}"}){_stencil_is_face_x(){
    {;_stencil_val(u.x,0,0,0); _stencil_val(u.x,0 -1,0,0);_stencil_val_a(uf.x,0,0,0);  }}end__stencil_is_face_x()
#line 199
_stencil_is_face_y(){
    {;_stencil_val(u.y,0,0,0); _stencil_val(u.y,0,0 -1,0);_stencil_val_a(uf.y,0,0,0);  }}end__stencil_is_face_y()
#line 199
_stencil_is_face_z(){
    {;_stencil_val(u.z,0,0,0); _stencil_val(u.z,0,0,0 -1);_stencil_val_a(uf.z,0,0,0);  }}end__stencil_is_face_z()}end_foreach_face_stencil();
  {
#line 199
foreach_face_generic(){is_face_x(){
    val(uf.x,0,0,0) = _const_fm.x*((val(u.x,0,0,0) + val(u.x,0 -1,0,0))/2.);}end_is_face_x()
#line 199
is_face_y(){
    val(uf.y,0,0,0) = _const_fm.y*((val(u.y,0,0,0) + val(u.y,0,0 -1,0))/2.);}end_is_face_y()
#line 199
is_face_z(){
    val(uf.z,0,0,0) = _const_fm.z*((val(u.z,0,0,0) + val(u.z,0,0,0 -1))/2.);}end_is_face_z()}end_foreach_face_generic();}}




  event ("properties");





  dtmax = DT;
  event ("stability");
}{end_tracing("init","/home/lisergey/basilisk/src/navier-stokes/centered.h",213);return 0;}end_tracing("init","/home/lisergey/basilisk/src/navier-stokes/centered.h",213);}








static int set_dtmax_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i++)!=0;*ip=i;*tp=t;return ret;}
#line 222 "/home/lisergey/basilisk/src/navier-stokes/centered.h"
      static int set_dtmax(const int i,const double t,Event *_ev){tracing("set_dtmax","/home/lisergey/basilisk/src/navier-stokes/centered.h",222); dtmax = DT;{end_tracing("set_dtmax","/home/lisergey/basilisk/src/navier-stokes/centered.h",222);return 0;}end_tracing("set_dtmax","/home/lisergey/basilisk/src/navier-stokes/centered.h",222);}

static int stability_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i++)!=0;*ip=i;*tp=t;return ret;}


#line 224
      static int stability(const int i,const double t,Event *_ev){tracing("stability","/home/lisergey/basilisk/src/navier-stokes/centered.h",224); {
  dt = dtnext (stokes ? dtmax : timestep (uf, dtmax));
}{end_tracing("stability","/home/lisergey/basilisk/src/navier-stokes/centered.h",226);return 0;}end_tracing("stability","/home/lisergey/basilisk/src/navier-stokes/centered.h",226);}







static int vof_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i++)!=0;*ip=i;*tp=t;return ret;}








#line 234
static int vof(const int i,const double t,Event *_ev){;return 0;}
static int tracer_advection_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i++)!=0;*ip=i;*tp=t;return ret;}

#line 235
static int tracer_advection(const int i,const double t,Event *_ev){;return 0;}
static int tracer_diffusion_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i++)!=0;*ip=i;*tp=t;return ret;}

#line 236
static int tracer_diffusion(const int i,const double t,Event *_ev){;return 0;}






static int properties_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i++)!=0;*ip=i;*tp=t;return ret;}







#line 243
static int properties(const int i,const double t,Event *_ev){;return 0;}
#line 255 "/home/lisergey/basilisk/src/navier-stokes/centered.h"
void prediction()
{
  vector du;
   {
    scalar s = new_scalar("s");
    du.x = s;
  } 
#line 258
{
    scalar s = new_scalar("s");
    du.y = s;
  } 
#line 258
{
    scalar s = new_scalar("s");
    du.z = s;
  }

  if (_attribute[u.x.i].gradient)
    {
    
#line 264
foreach_stencil(1,{(NonLocal[]){{"u","vector",(void *)&u,NULL,0},{"du","vector",(void *)&du,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n      \n// #line 265 \"/home/lisergey/basilisk/src/navier-stokes/centered.h\"\n{ {\n\n\n\n\n\n   val_out_(du.x,0,0,0) = u.x.gradient (val(u.x,-1,0,0), val(u.x,0,0,0), val(u.x,1,0,0))/Delta;\n      } \n// #line 265\n{\n\n\n\n\n\n   val_out_(du.y,0,0,0) = u.y.gradient (val(u.y,0,-1,0), val(u.y,0,0,0), val(u.y,0,1,0))/Delta;\n      } \n// #line 265\n{\n\n\n\n\n\n   val_out_(du.z,0,0,0) = u.z.gradient (val(u.z,0,0,-1), val(u.z,0,0,0), val(u.z,0,0,1))/Delta;\n      }}"})
      { {





_stencil_val(u.x,-1,0,0); _stencil_val(u.x,0,0,0); _stencil_val(u.x,1,0,0);





   
#line 271
_stencil_val_a(du.x,0,0,0);   
      } 
#line 265
{





_stencil_val(u.y,0,-1,0); _stencil_val(u.y,0,0,0); _stencil_val(u.y,0,1,0);





   
#line 271
_stencil_val_a(du.y,0,0,0);   
      } 
#line 265
{





_stencil_val(u.z,0,0,-1); _stencil_val(u.z,0,0,0); _stencil_val(u.z,0,0,1);





   
#line 271
_stencil_val_a(du.z,0,0,0);   
      }}end_foreach_stencil();{
#line 264
foreach()
      { {





   val(du.x,0,0,0) = _attribute[u.x.i].gradient (val(u.x,-1,0,0), val(u.x,0,0,0), val(u.x,1,0,0))/Delta;
      } 
#line 265
{





   val(du.y,0,0,0) = _attribute[u.y.i].gradient (val(u.y,0,-1,0), val(u.y,0,0,0), val(u.y,0,1,0))/Delta;
      } 
#line 265
{





   val(du.z,0,0,0) = _attribute[u.z.i].gradient (val(u.z,0,0,-1), val(u.z,0,0,0), val(u.z,0,0,1))/Delta;
      }}end_foreach();}}
  else
    {
    
#line 274
foreach_stencil(1,{(NonLocal[]){{"u","vector",(void *)&u,NULL,0},{"du","vector",(void *)&du,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n      \n// #line 275 \"/home/lisergey/basilisk/src/navier-stokes/centered.h\"\n{ {\n\n\n\n\n\n   val_out_(du.x,0,0,0) = (val(u.x,1,0,0) - val(u.x,-1,0,0))/(2.*Delta);\n    } \n// #line 275\n{\n\n\n\n\n\n   val_out_(du.y,0,0,0) = (val(u.y,0,1,0) - val(u.y,0,-1,0))/(2.*Delta);\n    } \n// #line 275\n{\n\n\n\n\n\n   val_out_(du.z,0,0,0) = (val(u.z,0,0,1) - val(u.z,0,0,-1))/(2.*Delta);\n    }}"})
      { {





_stencil_val(u.x,1,0,0); _stencil_val(u.x,-1,0,0);





   
#line 281
_stencil_val_a(du.x,0,0,0);   
    } 
#line 275
{





_stencil_val(u.y,0,1,0); _stencil_val(u.y,0,-1,0);





   
#line 281
_stencil_val_a(du.y,0,0,0);   
    } 
#line 275
{





_stencil_val(u.z,0,0,1); _stencil_val(u.z,0,0,-1);





   
#line 281
_stencil_val_a(du.z,0,0,0);   
    }}end_foreach_stencil();{
#line 274
foreach()
      { {





   val(du.x,0,0,0) = (val(u.x,1,0,0) - val(u.x,-1,0,0))/(2.*Delta);
    } 
#line 275
{





   val(du.y,0,0,0) = (val(u.y,0,1,0) - val(u.y,0,-1,0))/(2.*Delta);
    } 
#line 275
{





   val(du.z,0,0,0) = (val(u.z,0,0,1) - val(u.z,0,0,-1))/(2.*Delta);
    }}end_foreach();}}

  trash (((vector[]){uf,{{-1},{-1},{-1}}}));
  if(!is_constant(fm.x)){
  
#line 285
foreach_face_stencil(1,{(NonLocal[]){{"fm","vector",(void *)&fm,NULL,0},{"du","vector",(void *)&du,NULL,0},{"g","vector",(void *)&g,NULL,0},{"uf","vector",(void *)&uf,NULL,0},{"u","vector",(void *)&u,NULL,0},{"dt","double",(void *)&dt,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n// #line 285 \"/home/lisergey/basilisk/src/navier-stokes/centered.h\"\n{is_face_x(){ {\n    real un = dt*(val(u.x,0,0,0) + val(u.x,-1,0,0))/(2.*Delta), s = sign(un);\n    int i = -(s + 1.)/2.;\n    val_out_(uf.x,0,0,0) = val(u.x,i,0,0) + (val(g.x,0,0,0) + val(g.x,-1,0,0))*dt/4. + s*(1. - s*un)*val(du.x,i,0,0)*Delta/2.;\n\n    if (val(fm.y,i,0,0) && val(fm.y,i,1,0)) {\n      real fyy = val(u.y,i,0,0) < 0. ? val(u.x,i,1,0) - val(u.x,i,0,0) : val(u.x,i,0,0) - val(u.x,i,-1,0);\n      val_out_(uf.x,0,0,0) -= dt*val(u.y,i,0,0)*fyy/(2.*Delta);\n    }\n\n\n    if (val(fm.z,i,0,0) && val(fm.z,i,0,1)) {\n      real fzz = val(u.z,i,0,0) < 0. ? val(u.x,i,0,1) - val(u.x,i,0,0) : val(u.x,i,0,0) - val(u.x,i,0,-1);\n      val_out_(uf.x,0,0,0) -= dt*val(u.z,i,0,0)*fzz/(2.*Delta);\n    }\n\n    val_out_(uf.x,0,0,0) *= val(fm.x,0,0,0);\n  }}end_is_face_x()\n// #line 285\nis_face_y(){ {\n    real un = dt*(val(u.y,0,0,0) + val(u.y,0,-1,0))/(2.*Delta), s = sign(un);\n    int i = -(s + 1.)/2.;\n    val_out_(uf.y,0,0,0) = val(u.y,0,i,0) + (val(g.y,0,0,0) + val(g.y,0,-1,0))*dt/4. + s*(1. - s*un)*val(du.y,0,i,0)*Delta/2.;\n\n    if (val(fm.z,0,i,0) && val(fm.z,0,i,1)) {\n      real fyy = val(u.z,0,i,0) < 0. ? val(u.y,0,i,1) - val(u.y,0,i,0) : val(u.y,0,i,0) - val(u.y,0,i,-1);\n      val_out_(uf.y,0,0,0) -= dt*val(u.z,0,i,0)*fyy/(2.*Delta);\n    }\n\n\n    if (val(fm.x,0,i,0) && val(fm.x,1,i,0)) {\n      real fzz = val(u.x,0,i,0) < 0. ? val(u.y,1,i,0) - val(u.y,0,i,0) : val(u.y,0,i,0) - val(u.y,-1,i,0);\n      val_out_(uf.y,0,0,0) -= dt*val(u.x,0,i,0)*fzz/(2.*Delta);\n    }\n\n    val_out_(uf.y,0,0,0) *= val(fm.y,0,0,0);\n  }}end_is_face_y()\n// #line 285\nis_face_z(){ {\n    real un = dt*(val(u.z,0,0,0) + val(u.z,0,0,-1))/(2.*Delta), s = sign(un);\n    int i = -(s + 1.)/2.;\n    val_out_(uf.z,0,0,0) = val(u.z,0,0,i) + (val(g.z,0,0,0) + val(g.z,0,0,-1))*dt/4. + s*(1. - s*un)*val(du.z,0,0,i)*Delta/2.;\n\n    if (val(fm.x,0,0,i) && val(fm.x,1,0,i)) {\n      real fyy = val(u.x,0,0,i) < 0. ? val(u.z,1,0,i) - val(u.z,0,0,i) : val(u.z,0,0,i) - val(u.z,-1,0,i);\n      val_out_(uf.z,0,0,0) -= dt*val(u.x,0,0,i)*fyy/(2.*Delta);\n    }\n\n\n    if (val(fm.y,0,0,i) && val(fm.y,0,1,i)) {\n      real fzz = val(u.y,0,0,i) < 0. ? val(u.z,0,1,i) - val(u.z,0,0,i) : val(u.z,0,0,i) - val(u.z,0,-1,i);\n      val_out_(uf.z,0,0,0) -= dt*val(u.y,0,0,i)*fzz/(2.*Delta);\n    }\n\n    val_out_(uf.z,0,0,0) *= val(fm.z,0,0,0);\n  }}end_is_face_z()}"}){_stencil_is_face_x(){ {       
     _stencil_val(u.x,-1,0,0);_stencil_val(u.x,0,0,0);     
    
_stencil_val(u.x, o_stencil,0,0);_stencil_val(g.x,0,0,0); _stencil_val(g.x,-1,0,0);_stencil_val(du.x,o_stencil,0,0);
    
#line 288
_stencil_val_a(uf.x,0,0,0);

_stencil_val(fm.y,o_stencil,0,0); _stencil_val(fm.y,o_stencil,1,0); {        
       _stencil_val(u.x,o_stencil,-1,0);_stencil_val(u.x, o_stencil,0,0);_stencil_val(u.x, o_stencil,0,0); _stencil_val(u.x,o_stencil,1,0);_stencil_val(u.y, o_stencil,0,0);
_stencil_val(u.y,o_stencil,0,0);
      
#line 292
_stencil_val_r(uf.x,0,0,0);  
    }


_stencil_val(fm.z,o_stencil,0,0); _stencil_val(fm.z,o_stencil,0,1); {        
       _stencil_val(u.x,o_stencil,0,-1);_stencil_val(u.x, o_stencil,0,0);_stencil_val(u.x, o_stencil,0,0); _stencil_val(u.x,o_stencil,0,1);_stencil_val(u.z, o_stencil,0,0);
_stencil_val(u.z,o_stencil,0,0);
      
#line 298
_stencil_val_r(uf.x,0,0,0);  
    } 

_stencil_val(fm.x,0,0,0);        

      


      

    
#line 301
_stencil_val_r(uf.x,0,0,0); 
  }}end__stencil_is_face_x()
#line 285
_stencil_is_face_y(){ {       
     _stencil_val(u.y,0,-1,0);_stencil_val(u.y,0,0,0);     
    
_stencil_val(u.y,0, o_stencil,0);_stencil_val(g.y,0,0,0); _stencil_val(g.y,0,-1,0);_stencil_val(du.y,0,o_stencil,0);
    
#line 288
_stencil_val_a(uf.y,0,0,0);

_stencil_val(fm.z,0,o_stencil,0); _stencil_val(fm.z,0,o_stencil,1); {        
       _stencil_val(u.y,0,o_stencil,-1);_stencil_val(u.y,0, o_stencil,0);_stencil_val(u.y,0, o_stencil,0); _stencil_val(u.y,0,o_stencil,1);_stencil_val(u.z,0, o_stencil,0);
_stencil_val(u.z,0,o_stencil,0);
      
#line 292
_stencil_val_r(uf.y,0,0,0);  
    }


_stencil_val(fm.x,0,o_stencil,0); _stencil_val(fm.x,1,o_stencil,0); {        
       _stencil_val(u.y,-1,o_stencil,0);_stencil_val(u.y,0, o_stencil,0);_stencil_val(u.y,0, o_stencil,0); _stencil_val(u.y,1,o_stencil,0);_stencil_val(u.x,0, o_stencil,0);
_stencil_val(u.x,0,o_stencil,0);
      
#line 298
_stencil_val_r(uf.y,0,0,0);  
    } 

_stencil_val(fm.y,0,0,0);        

      


      

    
#line 301
_stencil_val_r(uf.y,0,0,0); 
  }}end__stencil_is_face_y()
#line 285
_stencil_is_face_z(){ {       
     _stencil_val(u.z,0,0,-1);_stencil_val(u.z,0,0,0);     
    
_stencil_val(u.z,0,0, o_stencil);_stencil_val(g.z,0,0,0); _stencil_val(g.z,0,0,-1);_stencil_val(du.z,0,0,o_stencil);
    
#line 288
_stencil_val_a(uf.z,0,0,0);

_stencil_val(fm.x,0,0,o_stencil); _stencil_val(fm.x,1,0,o_stencil); {        
       _stencil_val(u.z,-1,0,o_stencil);_stencil_val(u.z,0,0, o_stencil);_stencil_val(u.z,0,0, o_stencil); _stencil_val(u.z,1,0,o_stencil);_stencil_val(u.x,0,0, o_stencil);
_stencil_val(u.x,0,0,o_stencil);
      
#line 292
_stencil_val_r(uf.z,0,0,0);  
    }


_stencil_val(fm.y,0,0,o_stencil); _stencil_val(fm.y,0,1,o_stencil); {        
       _stencil_val(u.z,0,-1,o_stencil);_stencil_val(u.z,0,0, o_stencil);_stencil_val(u.z,0,0, o_stencil); _stencil_val(u.z,0,1,o_stencil);_stencil_val(u.y,0,0, o_stencil);
_stencil_val(u.y,0,0,o_stencil);
      
#line 298
_stencil_val_r(uf.z,0,0,0);  
    } 

_stencil_val(fm.z,0,0,0);        

      


      

    
#line 301
_stencil_val_r(uf.z,0,0,0); 
  }}end__stencil_is_face_z()}end_foreach_face_stencil();{
#line 285
foreach_face_generic(){is_face_x(){ {
    double un = dt*(val(u.x,0,0,0) + val(u.x,-1,0,0))/(2.*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    val(uf.x,0,0,0) = val(u.x,i,0,0) + (val(g.x,0,0,0) + val(g.x,-1,0,0))*dt/4. + s*(1. - s*un)*val(du.x,i,0,0)*Delta/2.;

    if (val(fm.y,i,0,0) && val(fm.y,i,1,0)) {
      double fyy = val(u.y,i,0,0) < 0. ? val(u.x,i,1,0) - val(u.x,i,0,0) : val(u.x,i,0,0) - val(u.x,i,-1,0);
      val(uf.x,0,0,0) -= dt*val(u.y,i,0,0)*fyy/(2.*Delta);
    }


    if (val(fm.z,i,0,0) && val(fm.z,i,0,1)) {
      double fzz = val(u.z,i,0,0) < 0. ? val(u.x,i,0,1) - val(u.x,i,0,0) : val(u.x,i,0,0) - val(u.x,i,0,-1);
      val(uf.x,0,0,0) -= dt*val(u.z,i,0,0)*fzz/(2.*Delta);
    }

    val(uf.x,0,0,0) *= val(fm.x,0,0,0);
  }}end_is_face_x()
#line 285
is_face_y(){ {
    double un = dt*(val(u.y,0,0,0) + val(u.y,0,-1,0))/(2.*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    val(uf.y,0,0,0) = val(u.y,0,i,0) + (val(g.y,0,0,0) + val(g.y,0,-1,0))*dt/4. + s*(1. - s*un)*val(du.y,0,i,0)*Delta/2.;

    if (val(fm.z,0,i,0) && val(fm.z,0,i,1)) {
      double fyy = val(u.z,0,i,0) < 0. ? val(u.y,0,i,1) - val(u.y,0,i,0) : val(u.y,0,i,0) - val(u.y,0,i,-1);
      val(uf.y,0,0,0) -= dt*val(u.z,0,i,0)*fyy/(2.*Delta);
    }


    if (val(fm.x,0,i,0) && val(fm.x,1,i,0)) {
      double fzz = val(u.x,0,i,0) < 0. ? val(u.y,1,i,0) - val(u.y,0,i,0) : val(u.y,0,i,0) - val(u.y,-1,i,0);
      val(uf.y,0,0,0) -= dt*val(u.x,0,i,0)*fzz/(2.*Delta);
    }

    val(uf.y,0,0,0) *= val(fm.y,0,0,0);
  }}end_is_face_y()
#line 285
is_face_z(){ {
    double un = dt*(val(u.z,0,0,0) + val(u.z,0,0,-1))/(2.*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    val(uf.z,0,0,0) = val(u.z,0,0,i) + (val(g.z,0,0,0) + val(g.z,0,0,-1))*dt/4. + s*(1. - s*un)*val(du.z,0,0,i)*Delta/2.;

    if (val(fm.x,0,0,i) && val(fm.x,1,0,i)) {
      double fyy = val(u.x,0,0,i) < 0. ? val(u.z,1,0,i) - val(u.z,0,0,i) : val(u.z,0,0,i) - val(u.z,-1,0,i);
      val(uf.z,0,0,0) -= dt*val(u.x,0,0,i)*fyy/(2.*Delta);
    }


    if (val(fm.y,0,0,i) && val(fm.y,0,1,i)) {
      double fzz = val(u.y,0,0,i) < 0. ? val(u.z,0,1,i) - val(u.z,0,0,i) : val(u.z,0,0,i) - val(u.z,0,-1,i);
      val(uf.z,0,0,0) -= dt*val(u.y,0,0,i)*fzz/(2.*Delta);
    }

    val(uf.z,0,0,0) *= val(fm.z,0,0,0);
  }}end_is_face_z()}end_foreach_face_generic();}}else {_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);
  
#line 285
foreach_face_stencil(1,{(NonLocal[]){{"_const_fm","_coord",(void *)&_const_fm,NULL,0},{"du","vector",(void *)&du,NULL,0},{"g","vector",(void *)&g,NULL,0},{"uf","vector",(void *)&uf,NULL,0},{"u","vector",(void *)&u,NULL,0},{"dt","double",(void *)&dt,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n// #line 285 \"/home/lisergey/basilisk/src/navier-stokes/centered.h\"\n{is_face_x(){ {\n    real un = dt*(val(u.x,0,0,0) + val(u.x,-1,0,0))/(2.*Delta), s = sign(un);\n    int i = -(s + 1.)/2.;\n    val_out_(uf.x,0,0,0) = val(u.x,i,0,0) + (val(g.x,0,0,0) + val(g.x,-1,0,0))*dt/4. + s*(1. - s*un)*val(du.x,i,0,0)*Delta/2.;\n\n    if (_const_fm.y && _const_fm.y) {\n      real fyy = val(u.y,i,0,0) < 0. ? val(u.x,i,1,0) - val(u.x,i,0,0) : val(u.x,i,0,0) - val(u.x,i,-1,0);\n      val_out_(uf.x,0,0,0) -= dt*val(u.y,i,0,0)*fyy/(2.*Delta);\n    }\n\n\n    if (_const_fm.z && _const_fm.z) {\n      real fzz = val(u.z,i,0,0) < 0. ? val(u.x,i,0,1) - val(u.x,i,0,0) : val(u.x,i,0,0) - val(u.x,i,0,-1);\n      val_out_(uf.x,0,0,0) -= dt*val(u.z,i,0,0)*fzz/(2.*Delta);\n    }\n\n    val_out_(uf.x,0,0,0) *= _const_fm.x;\n  }}end_is_face_x()\n// #line 285\nis_face_y(){ {\n    real un = dt*(val(u.y,0,0,0) + val(u.y,0,-1,0))/(2.*Delta), s = sign(un);\n    int i = -(s + 1.)/2.;\n    val_out_(uf.y,0,0,0) = val(u.y,0,i,0) + (val(g.y,0,0,0) + val(g.y,0,-1,0))*dt/4. + s*(1. - s*un)*val(du.y,0,i,0)*Delta/2.;\n\n    if (_const_fm.z && _const_fm.z) {\n      real fyy = val(u.z,0,i,0) < 0. ? val(u.y,0,i,1) - val(u.y,0,i,0) : val(u.y,0,i,0) - val(u.y,0,i,-1);\n      val_out_(uf.y,0,0,0) -= dt*val(u.z,0,i,0)*fyy/(2.*Delta);\n    }\n\n\n    if (_const_fm.x && _const_fm.x) {\n      real fzz = val(u.x,0,i,0) < 0. ? val(u.y,1,i,0) - val(u.y,0,i,0) : val(u.y,0,i,0) - val(u.y,-1,i,0);\n      val_out_(uf.y,0,0,0) -= dt*val(u.x,0,i,0)*fzz/(2.*Delta);\n    }\n\n    val_out_(uf.y,0,0,0) *= _const_fm.y;\n  }}end_is_face_y()\n// #line 285\nis_face_z(){ {\n    real un = dt*(val(u.z,0,0,0) + val(u.z,0,0,-1))/(2.*Delta), s = sign(un);\n    int i = -(s + 1.)/2.;\n    val_out_(uf.z,0,0,0) = val(u.z,0,0,i) + (val(g.z,0,0,0) + val(g.z,0,0,-1))*dt/4. + s*(1. - s*un)*val(du.z,0,0,i)*Delta/2.;\n\n    if (_const_fm.x && _const_fm.x) {\n      real fyy = val(u.x,0,0,i) < 0. ? val(u.z,1,0,i) - val(u.z,0,0,i) : val(u.z,0,0,i) - val(u.z,-1,0,i);\n      val_out_(uf.z,0,0,0) -= dt*val(u.x,0,0,i)*fyy/(2.*Delta);\n    }\n\n\n    if (_const_fm.y && _const_fm.y) {\n      real fzz = val(u.y,0,0,i) < 0. ? val(u.z,0,1,i) - val(u.z,0,0,i) : val(u.z,0,0,i) - val(u.z,0,-1,i);\n      val_out_(uf.z,0,0,0) -= dt*val(u.y,0,0,i)*fzz/(2.*Delta);\n    }\n\n    val_out_(uf.z,0,0,0) *= _const_fm.z;\n  }}end_is_face_z()}"}){_stencil_is_face_x(){ {       
     _stencil_val(u.x,-1,0,0);_stencil_val(u.x,0,0,0);     
    
_stencil_val(u.x, o_stencil,0,0);_stencil_val(g.x,0,0,0); _stencil_val(g.x,-1,0,0);_stencil_val(du.x,o_stencil,0,0);
    
#line 288
_stencil_val_a(uf.x,0,0,0);

;; {        
       _stencil_val(u.x,o_stencil,-1,0);_stencil_val(u.x, o_stencil,0,0);_stencil_val(u.x, o_stencil,0,0); _stencil_val(u.x,o_stencil,1,0);_stencil_val(u.y, o_stencil,0,0);
_stencil_val(u.y,o_stencil,0,0);
      
#line 292
_stencil_val_r(uf.x,0,0,0);  
    }


;; {        
       _stencil_val(u.x,o_stencil,0,-1);_stencil_val(u.x, o_stencil,0,0);_stencil_val(u.x, o_stencil,0,0); _stencil_val(u.x,o_stencil,0,1);_stencil_val(u.z, o_stencil,0,0);
_stencil_val(u.z,o_stencil,0,0);
      
#line 298
_stencil_val_r(uf.x,0,0,0);  
    }

;        

      


      

    
#line 301
_stencil_val_r(uf.x,0,0,0); 
  }}end__stencil_is_face_x()
#line 285
_stencil_is_face_y(){ {       
     _stencil_val(u.y,0,-1,0);_stencil_val(u.y,0,0,0);     
    
_stencil_val(u.y,0, o_stencil,0);_stencil_val(g.y,0,0,0); _stencil_val(g.y,0,-1,0);_stencil_val(du.y,0,o_stencil,0);
    
#line 288
_stencil_val_a(uf.y,0,0,0);

;; {        
       _stencil_val(u.y,0,o_stencil,-1);_stencil_val(u.y,0, o_stencil,0);_stencil_val(u.y,0, o_stencil,0); _stencil_val(u.y,0,o_stencil,1);_stencil_val(u.z,0, o_stencil,0);
_stencil_val(u.z,0,o_stencil,0);
      
#line 292
_stencil_val_r(uf.y,0,0,0);  
    }


;; {        
       _stencil_val(u.y,-1,o_stencil,0);_stencil_val(u.y,0, o_stencil,0);_stencil_val(u.y,0, o_stencil,0); _stencil_val(u.y,1,o_stencil,0);_stencil_val(u.x,0, o_stencil,0);
_stencil_val(u.x,0,o_stencil,0);
      
#line 298
_stencil_val_r(uf.y,0,0,0);  
    }

;        

      


      

    
#line 301
_stencil_val_r(uf.y,0,0,0); 
  }}end__stencil_is_face_y()
#line 285
_stencil_is_face_z(){ {       
     _stencil_val(u.z,0,0,-1);_stencil_val(u.z,0,0,0);     
    
_stencil_val(u.z,0,0, o_stencil);_stencil_val(g.z,0,0,0); _stencil_val(g.z,0,0,-1);_stencil_val(du.z,0,0,o_stencil);
    
#line 288
_stencil_val_a(uf.z,0,0,0);

;; {        
       _stencil_val(u.z,-1,0,o_stencil);_stencil_val(u.z,0,0, o_stencil);_stencil_val(u.z,0,0, o_stencil); _stencil_val(u.z,1,0,o_stencil);_stencil_val(u.x,0,0, o_stencil);
_stencil_val(u.x,0,0,o_stencil);
      
#line 292
_stencil_val_r(uf.z,0,0,0);  
    }


;; {        
       _stencil_val(u.z,0,-1,o_stencil);_stencil_val(u.z,0,0, o_stencil);_stencil_val(u.z,0,0, o_stencil); _stencil_val(u.z,0,1,o_stencil);_stencil_val(u.y,0,0, o_stencil);
_stencil_val(u.y,0,0,o_stencil);
      
#line 298
_stencil_val_r(uf.z,0,0,0);  
    }

;        

      


      

    
#line 301
_stencil_val_r(uf.z,0,0,0); 
  }}end__stencil_is_face_z()}end_foreach_face_stencil();
  {
#line 285
foreach_face_generic(){is_face_x(){ {
    double un = dt*(val(u.x,0,0,0) + val(u.x,-1,0,0))/(2.*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    val(uf.x,0,0,0) = val(u.x,i,0,0) + (val(g.x,0,0,0) + val(g.x,-1,0,0))*dt/4. + s*(1. - s*un)*val(du.x,i,0,0)*Delta/2.;

    if (_const_fm.y && _const_fm.y) {
      double fyy = val(u.y,i,0,0) < 0. ? val(u.x,i,1,0) - val(u.x,i,0,0) : val(u.x,i,0,0) - val(u.x,i,-1,0);
      val(uf.x,0,0,0) -= dt*val(u.y,i,0,0)*fyy/(2.*Delta);
    }


    if (_const_fm.z && _const_fm.z) {
      double fzz = val(u.z,i,0,0) < 0. ? val(u.x,i,0,1) - val(u.x,i,0,0) : val(u.x,i,0,0) - val(u.x,i,0,-1);
      val(uf.x,0,0,0) -= dt*val(u.z,i,0,0)*fzz/(2.*Delta);
    }

    val(uf.x,0,0,0) *= _const_fm.x;
  }}end_is_face_x()
#line 285
is_face_y(){ {
    double un = dt*(val(u.y,0,0,0) + val(u.y,0,-1,0))/(2.*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    val(uf.y,0,0,0) = val(u.y,0,i,0) + (val(g.y,0,0,0) + val(g.y,0,-1,0))*dt/4. + s*(1. - s*un)*val(du.y,0,i,0)*Delta/2.;

    if (_const_fm.z && _const_fm.z) {
      double fyy = val(u.z,0,i,0) < 0. ? val(u.y,0,i,1) - val(u.y,0,i,0) : val(u.y,0,i,0) - val(u.y,0,i,-1);
      val(uf.y,0,0,0) -= dt*val(u.z,0,i,0)*fyy/(2.*Delta);
    }


    if (_const_fm.x && _const_fm.x) {
      double fzz = val(u.x,0,i,0) < 0. ? val(u.y,1,i,0) - val(u.y,0,i,0) : val(u.y,0,i,0) - val(u.y,-1,i,0);
      val(uf.y,0,0,0) -= dt*val(u.x,0,i,0)*fzz/(2.*Delta);
    }

    val(uf.y,0,0,0) *= _const_fm.y;
  }}end_is_face_y()
#line 285
is_face_z(){ {
    double un = dt*(val(u.z,0,0,0) + val(u.z,0,0,-1))/(2.*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    val(uf.z,0,0,0) = val(u.z,0,0,i) + (val(g.z,0,0,0) + val(g.z,0,0,-1))*dt/4. + s*(1. - s*un)*val(du.z,0,0,i)*Delta/2.;

    if (_const_fm.x && _const_fm.x) {
      double fyy = val(u.x,0,0,i) < 0. ? val(u.z,1,0,i) - val(u.z,0,0,i) : val(u.z,0,0,i) - val(u.z,-1,0,i);
      val(uf.z,0,0,0) -= dt*val(u.x,0,0,i)*fyy/(2.*Delta);
    }


    if (_const_fm.y && _const_fm.y) {
      double fzz = val(u.y,0,0,i) < 0. ? val(u.z,0,1,i) - val(u.z,0,0,i) : val(u.z,0,0,i) - val(u.z,0,-1,i);
      val(uf.z,0,0,0) -= dt*val(u.y,0,0,i)*fzz/(2.*Delta);
    }

    val(uf.z,0,0,0) *= _const_fm.z;
  }}end_is_face_z()}end_foreach_face_generic();}}

  delete ((scalar *)((vector[]){du,{{-1},{-1},{-1}}}));
}
#line 316
static int advection_term_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i++)!=0;*ip=i;*tp=t;return ret;}
#line 316 "/home/lisergey/basilisk/src/navier-stokes/centered.h"
      static int advection_term(const int i,const double t,Event *_ev){tracing("advection_term","/home/lisergey/basilisk/src/navier-stokes/centered.h",316);
{
  if (!stokes) {
    prediction();
    mgpf = project (uf, pf, alpha, dt/2., mgpf.nrelax);
    advection ((scalar *)((vector[]){u,{{-1},{-1},{-1}}}), uf, dt, (scalar *)((vector[]){g,{{-1},{-1},{-1}}}));
  }
}{end_tracing("advection_term","/home/lisergey/basilisk/src/navier-stokes/centered.h",323);return 0;}end_tracing("advection_term","/home/lisergey/basilisk/src/navier-stokes/centered.h",323);}







static void correction (double dt)
{
  foreach_stencil(1,{(NonLocal[]){{"g","vector",(void *)&g,NULL,0},{"dt","double",(void *)&dt,NULL,0},{"u","vector",(void *)&u,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n    \n// #line 334 \"/home/lisergey/basilisk/src/navier-stokes/centered.h\"\n{\n      val_out_(u.x,0,0,0) += dt*val(g.x,0,0,0);\n      \n// #line 335\nval_out_(u.y,0,0,0) += dt*val(g.y,0,0,0);\n      \n// #line 335\nval_out_(u.z,0,0,0) += dt*val(g.z,0,0,0);}"})
    {
      {_stencil_val(g.x,0,0,0);_stencil_val_r(u.x,0,0,0);  }
      
#line 335
{_stencil_val(g.y,0,0,0);_stencil_val_r(u.y,0,0,0);  }
      
#line 335
{_stencil_val(g.z,0,0,0);_stencil_val_r(u.z,0,0,0);  }}end_foreach_stencil();
  {
#line 333
foreach()
    {
      val(u.x,0,0,0) += dt*val(g.x,0,0,0);
      
#line 335
val(u.y,0,0,0) += dt*val(g.y,0,0,0);
      
#line 335
val(u.z,0,0,0) += dt*val(g.z,0,0,0);}end_foreach();}
}








static int viscous_term_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i++)!=0;*ip=i;*tp=t;return ret;}
#line 345 "/home/lisergey/basilisk/src/navier-stokes/centered.h"
      static int viscous_term(const int i,const double t,Event *_ev){tracing("viscous_term","/home/lisergey/basilisk/src/navier-stokes/centered.h",345);
{
  if (constant(mu.x) != 0.) {
    correction (dt);
    mgu = viscosity (u, mu, rho, dt, mgu.nrelax
#line 290 "/home/lisergey/basilisk/src/viscosity.h"
, NULL
#line 349 "/home/lisergey/basilisk/src/navier-stokes/centered.h"
);
    correction (-dt);
  }




  if (!is_constant(a.x)) {
    vector af = a;
    trash (((vector[]){af,{{-1},{-1},{-1}}}));
    foreach_face_stencil(1,{(NonLocal[]){{"af","vector",(void *)&af,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n// #line 359 \"/home/lisergey/basilisk/src/navier-stokes/centered.h\"\n{is_face_x(){\n      val_out_(af.x,0,0,0) = 0.;}end_is_face_x()\n// #line 359\nis_face_y(){\n      val_out_(af.y,0,0,0) = 0.;}end_is_face_y()\n// #line 359\nis_face_z(){\n      val_out_(af.z,0,0,0) = 0.;}end_is_face_z()}"}){_stencil_is_face_x(){
      {_stencil_val_a(af.x,0,0,0);  }}end__stencil_is_face_x()
#line 359
_stencil_is_face_y(){
      {_stencil_val_a(af.y,0,0,0);  }}end__stencil_is_face_y()
#line 359
_stencil_is_face_z(){
      {_stencil_val_a(af.z,0,0,0);  }}end__stencil_is_face_z()}end_foreach_face_stencil();
    {
#line 359
foreach_face_generic(){is_face_x(){
      val(af.x,0,0,0) = 0.;}end_is_face_x()
#line 359
is_face_y(){
      val(af.y,0,0,0) = 0.;}end_is_face_y()
#line 359
is_face_z(){
      val(af.z,0,0,0) = 0.;}end_is_face_z()}end_foreach_face_generic();}
  }
}{end_tracing("viscous_term","/home/lisergey/basilisk/src/navier-stokes/centered.h",362);return 0;}end_tracing("viscous_term","/home/lisergey/basilisk/src/navier-stokes/centered.h",362);}
#line 381
static int acceleration_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i++)!=0;*ip=i;*tp=t;return ret;}
#line 381 "/home/lisergey/basilisk/src/navier-stokes/centered.h"
      static int acceleration(const int i,const double t,Event *_ev){tracing("acceleration","/home/lisergey/basilisk/src/navier-stokes/centered.h",381);
{
  trash (((vector[]){uf,{{-1},{-1},{-1}}}));
  if(!is_constant(fm.x) && !is_constant(a.x)){
  
#line 384
foreach_face_stencil(1,{(NonLocal[]){{"a","vector",(void *)&a,NULL,0},{"dt","double",(void *)&dt,NULL,0},{"u","vector",(void *)&u,NULL,0},{"fm","vector",(void *)&fm,NULL,0},{"uf","vector",(void *)&uf,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n// #line 384 \"/home/lisergey/basilisk/src/navier-stokes/centered.h\"\n{is_face_x(){\n    val_out_(uf.x,0,0,0) = val(fm.x,0,0,0)*(((val(u.x,0,0,0) + val(u.x,0 -1,0,0))/2.) + dt*val(a.x,0,0,0));}end_is_face_x()\n// #line 384\nis_face_y(){\n    val_out_(uf.y,0,0,0) = val(fm.y,0,0,0)*(((val(u.y,0,0,0) + val(u.y,0,0 -1,0))/2.) + dt*val(a.y,0,0,0));}end_is_face_y()\n// #line 384\nis_face_z(){\n    val_out_(uf.z,0,0,0) = val(fm.z,0,0,0)*(((val(u.z,0,0,0) + val(u.z,0,0,0 -1))/2.) + dt*val(a.z,0,0,0));}end_is_face_z()}"}){_stencil_is_face_x(){
    { _stencil_val(fm.x,0,0,0);_stencil_val(u.x,0,0,0); _stencil_val(u.x,0 -1,0,0);_stencil_val(a.x,0,0,0);_stencil_val_a(uf.x,0,0,0);    }}end__stencil_is_face_x()
#line 384
_stencil_is_face_y(){
    { _stencil_val(fm.y,0,0,0);_stencil_val(u.y,0,0,0); _stencil_val(u.y,0,0 -1,0);_stencil_val(a.y,0,0,0);_stencil_val_a(uf.y,0,0,0);    }}end__stencil_is_face_y()
#line 384
_stencil_is_face_z(){
    { _stencil_val(fm.z,0,0,0);_stencil_val(u.z,0,0,0); _stencil_val(u.z,0,0,0 -1);_stencil_val(a.z,0,0,0);_stencil_val_a(uf.z,0,0,0);    }}end__stencil_is_face_z()}end_foreach_face_stencil();{
#line 384
foreach_face_generic(){is_face_x(){
    val(uf.x,0,0,0) = val(fm.x,0,0,0)*(((val(u.x,0,0,0) + val(u.x,0 -1,0,0))/2.) + dt*val(a.x,0,0,0));}end_is_face_x()
#line 384
is_face_y(){
    val(uf.y,0,0,0) = val(fm.y,0,0,0)*(((val(u.y,0,0,0) + val(u.y,0,0 -1,0))/2.) + dt*val(a.y,0,0,0));}end_is_face_y()
#line 384
is_face_z(){
    val(uf.z,0,0,0) = val(fm.z,0,0,0)*(((val(u.z,0,0,0) + val(u.z,0,0,0 -1))/2.) + dt*val(a.z,0,0,0));}end_is_face_z()}end_foreach_face_generic();}}else if(is_constant(fm.x) && !is_constant(a.x)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);
  
#line 384
foreach_face_stencil(1,{(NonLocal[]){{"a","vector",(void *)&a,NULL,0},{"dt","double",(void *)&dt,NULL,0},{"u","vector",(void *)&u,NULL,0},{"_const_fm","_coord",(void *)&_const_fm,NULL,0},{"uf","vector",(void *)&uf,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n// #line 384 \"/home/lisergey/basilisk/src/navier-stokes/centered.h\"\n{is_face_x(){\n    val_out_(uf.x,0,0,0) = _const_fm.x*(((val(u.x,0,0,0) + val(u.x,0 -1,0,0))/2.) + dt*val(a.x,0,0,0));}end_is_face_x()\n// #line 384\nis_face_y(){\n    val_out_(uf.y,0,0,0) = _const_fm.y*(((val(u.y,0,0,0) + val(u.y,0,0 -1,0))/2.) + dt*val(a.y,0,0,0));}end_is_face_y()\n// #line 384\nis_face_z(){\n    val_out_(uf.z,0,0,0) = _const_fm.z*(((val(u.z,0,0,0) + val(u.z,0,0,0 -1))/2.) + dt*val(a.z,0,0,0));}end_is_face_z()}"}){_stencil_is_face_x(){
    {;_stencil_val(u.x,0,0,0); _stencil_val(u.x,0 -1,0,0);_stencil_val(a.x,0,0,0);_stencil_val_a(uf.x,0,0,0);    }}end__stencil_is_face_x()
#line 384
_stencil_is_face_y(){
    {;_stencil_val(u.y,0,0,0); _stencil_val(u.y,0,0 -1,0);_stencil_val(a.y,0,0,0);_stencil_val_a(uf.y,0,0,0);    }}end__stencil_is_face_y()
#line 384
_stencil_is_face_z(){
    {;_stencil_val(u.z,0,0,0); _stencil_val(u.z,0,0,0 -1);_stencil_val(a.z,0,0,0);_stencil_val_a(uf.z,0,0,0);    }}end__stencil_is_face_z()}end_foreach_face_stencil();
  {
#line 384
foreach_face_generic(){is_face_x(){
    val(uf.x,0,0,0) = _const_fm.x*(((val(u.x,0,0,0) + val(u.x,0 -1,0,0))/2.) + dt*val(a.x,0,0,0));}end_is_face_x()
#line 384
is_face_y(){
    val(uf.y,0,0,0) = _const_fm.y*(((val(u.y,0,0,0) + val(u.y,0,0 -1,0))/2.) + dt*val(a.y,0,0,0));}end_is_face_y()
#line 384
is_face_z(){
    val(uf.z,0,0,0) = _const_fm.z*(((val(u.z,0,0,0) + val(u.z,0,0,0 -1))/2.) + dt*val(a.z,0,0,0));}end_is_face_z()}end_foreach_face_generic();}}else if(!is_constant(fm.x) && is_constant(a.x)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX],_constant[a.z.i-_NVARMAX]};NOT_UNUSED(_const_a);
  
#line 384
foreach_face_stencil(1,{(NonLocal[]){{"_const_a","_coord",(void *)&_const_a,NULL,0},{"dt","double",(void *)&dt,NULL,0},{"u","vector",(void *)&u,NULL,0},{"fm","vector",(void *)&fm,NULL,0},{"uf","vector",(void *)&uf,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n// #line 384 \"/home/lisergey/basilisk/src/navier-stokes/centered.h\"\n{is_face_x(){\n    val_out_(uf.x,0,0,0) = val(fm.x,0,0,0)*(((val(u.x,0,0,0) + val(u.x,0 -1,0,0))/2.) + dt*_const_a.x);}end_is_face_x()\n// #line 384\nis_face_y(){\n    val_out_(uf.y,0,0,0) = val(fm.y,0,0,0)*(((val(u.y,0,0,0) + val(u.y,0,0 -1,0))/2.) + dt*_const_a.y);}end_is_face_y()\n// #line 384\nis_face_z(){\n    val_out_(uf.z,0,0,0) = val(fm.z,0,0,0)*(((val(u.z,0,0,0) + val(u.z,0,0,0 -1))/2.) + dt*_const_a.z);}end_is_face_z()}"}){_stencil_is_face_x(){
    { _stencil_val(fm.x,0,0,0);_stencil_val(u.x,0,0,0); _stencil_val(u.x,0 -1,0,0);;_stencil_val_a(uf.x,0,0,0);    }}end__stencil_is_face_x()
#line 384
_stencil_is_face_y(){
    { _stencil_val(fm.y,0,0,0);_stencil_val(u.y,0,0,0); _stencil_val(u.y,0,0 -1,0);;_stencil_val_a(uf.y,0,0,0);    }}end__stencil_is_face_y()
#line 384
_stencil_is_face_z(){
    { _stencil_val(fm.z,0,0,0);_stencil_val(u.z,0,0,0); _stencil_val(u.z,0,0,0 -1);;_stencil_val_a(uf.z,0,0,0);    }}end__stencil_is_face_z()}end_foreach_face_stencil();
  {
#line 384
foreach_face_generic(){is_face_x(){
    val(uf.x,0,0,0) = val(fm.x,0,0,0)*(((val(u.x,0,0,0) + val(u.x,0 -1,0,0))/2.) + dt*_const_a.x);}end_is_face_x()
#line 384
is_face_y(){
    val(uf.y,0,0,0) = val(fm.y,0,0,0)*(((val(u.y,0,0,0) + val(u.y,0,0 -1,0))/2.) + dt*_const_a.y);}end_is_face_y()
#line 384
is_face_z(){
    val(uf.z,0,0,0) = val(fm.z,0,0,0)*(((val(u.z,0,0,0) + val(u.z,0,0,0 -1))/2.) + dt*_const_a.z);}end_is_face_z()}end_foreach_face_generic();}}else {_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX],_constant[a.z.i-_NVARMAX]};NOT_UNUSED(_const_a);
  
#line 384
foreach_face_stencil(1,{(NonLocal[]){{"_const_a","_coord",(void *)&_const_a,NULL,0},{"dt","double",(void *)&dt,NULL,0},{"u","vector",(void *)&u,NULL,0},{"_const_fm","_coord",(void *)&_const_fm,NULL,0},{"uf","vector",(void *)&uf,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n// #line 384 \"/home/lisergey/basilisk/src/navier-stokes/centered.h\"\n{is_face_x(){\n    val_out_(uf.x,0,0,0) = _const_fm.x*(((val(u.x,0,0,0) + val(u.x,0 -1,0,0))/2.) + dt*_const_a.x);}end_is_face_x()\n// #line 384\nis_face_y(){\n    val_out_(uf.y,0,0,0) = _const_fm.y*(((val(u.y,0,0,0) + val(u.y,0,0 -1,0))/2.) + dt*_const_a.y);}end_is_face_y()\n// #line 384\nis_face_z(){\n    val_out_(uf.z,0,0,0) = _const_fm.z*(((val(u.z,0,0,0) + val(u.z,0,0,0 -1))/2.) + dt*_const_a.z);}end_is_face_z()}"}){_stencil_is_face_x(){
    {;_stencil_val(u.x,0,0,0); _stencil_val(u.x,0 -1,0,0);;_stencil_val_a(uf.x,0,0,0);    }}end__stencil_is_face_x()
#line 384
_stencil_is_face_y(){
    {;_stencil_val(u.y,0,0,0); _stencil_val(u.y,0,0 -1,0);;_stencil_val_a(uf.y,0,0,0);    }}end__stencil_is_face_y()
#line 384
_stencil_is_face_z(){
    {;_stencil_val(u.z,0,0,0); _stencil_val(u.z,0,0,0 -1);;_stencil_val_a(uf.z,0,0,0);    }}end__stencil_is_face_z()}end_foreach_face_stencil();
  {
#line 384
foreach_face_generic(){is_face_x(){
    val(uf.x,0,0,0) = _const_fm.x*(((val(u.x,0,0,0) + val(u.x,0 -1,0,0))/2.) + dt*_const_a.x);}end_is_face_x()
#line 384
is_face_y(){
    val(uf.y,0,0,0) = _const_fm.y*(((val(u.y,0,0,0) + val(u.y,0,0 -1,0))/2.) + dt*_const_a.y);}end_is_face_y()
#line 384
is_face_z(){
    val(uf.z,0,0,0) = _const_fm.z*(((val(u.z,0,0,0) + val(u.z,0,0,0 -1))/2.) + dt*_const_a.z);}end_is_face_z()}end_foreach_face_generic();}}
}{end_tracing("acceleration","/home/lisergey/basilisk/src/navier-stokes/centered.h",386);return 0;}end_tracing("acceleration","/home/lisergey/basilisk/src/navier-stokes/centered.h",386);}
#line 395 "/home/lisergey/basilisk/src/navier-stokes/centered.h"
void centered_gradient (scalar p, vector g)
{





  vector  gf=new_face_vector("gf");
  if(!is_constant(fm.x) && !is_constant(a.x) && !is_constant(alpha.x)){
  
#line 403
foreach_face_stencil(1,{(NonLocal[]){{"p","scalar",(void *)&p,NULL,0},{"alpha","vector",(void *)&alpha,NULL,0},{"a","vector",(void *)&a,NULL,0},{"fm","vector",(void *)&fm,NULL,0},{"gf","vector",(void *)&gf,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n// #line 403 \"/home/lisergey/basilisk/src/navier-stokes/centered.h\"\n{is_face_x(){\n    val_out_(gf.x,0,0,0) = val(fm.x,0,0,0)*val(a.x,0,0,0) - val(alpha.x,0,0,0)*(val(p,0,0,0) - val(p,-1,0,0))/Delta;}end_is_face_x()\n// #line 403\nis_face_y(){\n    val_out_(gf.y,0,0,0) = val(fm.y,0,0,0)*val(a.y,0,0,0) - val(alpha.y,0,0,0)*(val(p,0,0,0) - val(p,0,-1,0))/Delta;}end_is_face_y()\n// #line 403\nis_face_z(){\n    val_out_(gf.z,0,0,0) = val(fm.z,0,0,0)*val(a.z,0,0,0) - val(alpha.z,0,0,0)*(val(p,0,0,0) - val(p,0,0,-1))/Delta;}end_is_face_z()}"}){_stencil_is_face_x(){
    { _stencil_val(fm.x,0,0,0);_stencil_val(a.x,0,0,0); _stencil_val(alpha.x,0,0,0);_stencil_val(p,0,0,0); _stencil_val(p,-1,0,0);_stencil_val_a(gf.x,0,0,0);   }}end__stencil_is_face_x()
#line 403
_stencil_is_face_y(){
    { _stencil_val(fm.y,0,0,0);_stencil_val(a.y,0,0,0); _stencil_val(alpha.y,0,0,0);_stencil_val(p,0,0,0); _stencil_val(p,0,-1,0);_stencil_val_a(gf.y,0,0,0);   }}end__stencil_is_face_y()
#line 403
_stencil_is_face_z(){
    { _stencil_val(fm.z,0,0,0);_stencil_val(a.z,0,0,0); _stencil_val(alpha.z,0,0,0);_stencil_val(p,0,0,0); _stencil_val(p,0,0,-1);_stencil_val_a(gf.z,0,0,0);   }}end__stencil_is_face_z()}end_foreach_face_stencil();{
#line 403
foreach_face_generic(){is_face_x(){
    val(gf.x,0,0,0) = val(fm.x,0,0,0)*val(a.x,0,0,0) - val(alpha.x,0,0,0)*(val(p,0,0,0) - val(p,-1,0,0))/Delta;}end_is_face_x()
#line 403
is_face_y(){
    val(gf.y,0,0,0) = val(fm.y,0,0,0)*val(a.y,0,0,0) - val(alpha.y,0,0,0)*(val(p,0,0,0) - val(p,0,-1,0))/Delta;}end_is_face_y()
#line 403
is_face_z(){
    val(gf.z,0,0,0) = val(fm.z,0,0,0)*val(a.z,0,0,0) - val(alpha.z,0,0,0)*(val(p,0,0,0) - val(p,0,0,-1))/Delta;}end_is_face_z()}end_foreach_face_generic();}}else if(is_constant(fm.x) && !is_constant(a.x) && !is_constant(alpha.x)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);
  
#line 403
foreach_face_stencil(1,{(NonLocal[]){{"p","scalar",(void *)&p,NULL,0},{"alpha","vector",(void *)&alpha,NULL,0},{"a","vector",(void *)&a,NULL,0},{"_const_fm","_coord",(void *)&_const_fm,NULL,0},{"gf","vector",(void *)&gf,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n// #line 403 \"/home/lisergey/basilisk/src/navier-stokes/centered.h\"\n{is_face_x(){\n    val_out_(gf.x,0,0,0) = _const_fm.x*val(a.x,0,0,0) - val(alpha.x,0,0,0)*(val(p,0,0,0) - val(p,-1,0,0))/Delta;}end_is_face_x()\n// #line 403\nis_face_y(){\n    val_out_(gf.y,0,0,0) = _const_fm.y*val(a.y,0,0,0) - val(alpha.y,0,0,0)*(val(p,0,0,0) - val(p,0,-1,0))/Delta;}end_is_face_y()\n// #line 403\nis_face_z(){\n    val_out_(gf.z,0,0,0) = _const_fm.z*val(a.z,0,0,0) - val(alpha.z,0,0,0)*(val(p,0,0,0) - val(p,0,0,-1))/Delta;}end_is_face_z()}"}){_stencil_is_face_x(){
    {;_stencil_val(a.x,0,0,0); _stencil_val(alpha.x,0,0,0);_stencil_val(p,0,0,0); _stencil_val(p,-1,0,0);_stencil_val_a(gf.x,0,0,0);   }}end__stencil_is_face_x()
#line 403
_stencil_is_face_y(){
    {;_stencil_val(a.y,0,0,0); _stencil_val(alpha.y,0,0,0);_stencil_val(p,0,0,0); _stencil_val(p,0,-1,0);_stencil_val_a(gf.y,0,0,0);   }}end__stencil_is_face_y()
#line 403
_stencil_is_face_z(){
    {;_stencil_val(a.z,0,0,0); _stencil_val(alpha.z,0,0,0);_stencil_val(p,0,0,0); _stencil_val(p,0,0,-1);_stencil_val_a(gf.z,0,0,0);   }}end__stencil_is_face_z()}end_foreach_face_stencil();
  {
#line 403
foreach_face_generic(){is_face_x(){
    val(gf.x,0,0,0) = _const_fm.x*val(a.x,0,0,0) - val(alpha.x,0,0,0)*(val(p,0,0,0) - val(p,-1,0,0))/Delta;}end_is_face_x()
#line 403
is_face_y(){
    val(gf.y,0,0,0) = _const_fm.y*val(a.y,0,0,0) - val(alpha.y,0,0,0)*(val(p,0,0,0) - val(p,0,-1,0))/Delta;}end_is_face_y()
#line 403
is_face_z(){
    val(gf.z,0,0,0) = _const_fm.z*val(a.z,0,0,0) - val(alpha.z,0,0,0)*(val(p,0,0,0) - val(p,0,0,-1))/Delta;}end_is_face_z()}end_foreach_face_generic();}}else if(!is_constant(fm.x) && is_constant(a.x) && !is_constant(alpha.x)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX],_constant[a.z.i-_NVARMAX]};NOT_UNUSED(_const_a);
  
#line 403
foreach_face_stencil(1,{(NonLocal[]){{"p","scalar",(void *)&p,NULL,0},{"alpha","vector",(void *)&alpha,NULL,0},{"_const_a","_coord",(void *)&_const_a,NULL,0},{"fm","vector",(void *)&fm,NULL,0},{"gf","vector",(void *)&gf,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n// #line 403 \"/home/lisergey/basilisk/src/navier-stokes/centered.h\"\n{is_face_x(){\n    val_out_(gf.x,0,0,0) = val(fm.x,0,0,0)*_const_a.x - val(alpha.x,0,0,0)*(val(p,0,0,0) - val(p,-1,0,0))/Delta;}end_is_face_x()\n// #line 403\nis_face_y(){\n    val_out_(gf.y,0,0,0) = val(fm.y,0,0,0)*_const_a.y - val(alpha.y,0,0,0)*(val(p,0,0,0) - val(p,0,-1,0))/Delta;}end_is_face_y()\n// #line 403\nis_face_z(){\n    val_out_(gf.z,0,0,0) = val(fm.z,0,0,0)*_const_a.z - val(alpha.z,0,0,0)*(val(p,0,0,0) - val(p,0,0,-1))/Delta;}end_is_face_z()}"}){_stencil_is_face_x(){
    { _stencil_val(fm.x,0,0,0);; _stencil_val(alpha.x,0,0,0);_stencil_val(p,0,0,0); _stencil_val(p,-1,0,0);_stencil_val_a(gf.x,0,0,0);   }}end__stencil_is_face_x()
#line 403
_stencil_is_face_y(){
    { _stencil_val(fm.y,0,0,0);; _stencil_val(alpha.y,0,0,0);_stencil_val(p,0,0,0); _stencil_val(p,0,-1,0);_stencil_val_a(gf.y,0,0,0);   }}end__stencil_is_face_y()
#line 403
_stencil_is_face_z(){
    { _stencil_val(fm.z,0,0,0);; _stencil_val(alpha.z,0,0,0);_stencil_val(p,0,0,0); _stencil_val(p,0,0,-1);_stencil_val_a(gf.z,0,0,0);   }}end__stencil_is_face_z()}end_foreach_face_stencil();
  {
#line 403
foreach_face_generic(){is_face_x(){
    val(gf.x,0,0,0) = val(fm.x,0,0,0)*_const_a.x - val(alpha.x,0,0,0)*(val(p,0,0,0) - val(p,-1,0,0))/Delta;}end_is_face_x()
#line 403
is_face_y(){
    val(gf.y,0,0,0) = val(fm.y,0,0,0)*_const_a.y - val(alpha.y,0,0,0)*(val(p,0,0,0) - val(p,0,-1,0))/Delta;}end_is_face_y()
#line 403
is_face_z(){
    val(gf.z,0,0,0) = val(fm.z,0,0,0)*_const_a.z - val(alpha.z,0,0,0)*(val(p,0,0,0) - val(p,0,0,-1))/Delta;}end_is_face_z()}end_foreach_face_generic();}}else if(is_constant(fm.x) && is_constant(a.x) && !is_constant(alpha.x)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX],_constant[a.z.i-_NVARMAX]};NOT_UNUSED(_const_a);
  
#line 403
foreach_face_stencil(1,{(NonLocal[]){{"p","scalar",(void *)&p,NULL,0},{"alpha","vector",(void *)&alpha,NULL,0},{"_const_a","_coord",(void *)&_const_a,NULL,0},{"_const_fm","_coord",(void *)&_const_fm,NULL,0},{"gf","vector",(void *)&gf,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n// #line 403 \"/home/lisergey/basilisk/src/navier-stokes/centered.h\"\n{is_face_x(){\n    val_out_(gf.x,0,0,0) = _const_fm.x*_const_a.x - val(alpha.x,0,0,0)*(val(p,0,0,0) - val(p,-1,0,0))/Delta;}end_is_face_x()\n// #line 403\nis_face_y(){\n    val_out_(gf.y,0,0,0) = _const_fm.y*_const_a.y - val(alpha.y,0,0,0)*(val(p,0,0,0) - val(p,0,-1,0))/Delta;}end_is_face_y()\n// #line 403\nis_face_z(){\n    val_out_(gf.z,0,0,0) = _const_fm.z*_const_a.z - val(alpha.z,0,0,0)*(val(p,0,0,0) - val(p,0,0,-1))/Delta;}end_is_face_z()}"}){_stencil_is_face_x(){
    {;; _stencil_val(alpha.x,0,0,0);_stencil_val(p,0,0,0); _stencil_val(p,-1,0,0);_stencil_val_a(gf.x,0,0,0);   }}end__stencil_is_face_x()
#line 403
_stencil_is_face_y(){
    {;; _stencil_val(alpha.y,0,0,0);_stencil_val(p,0,0,0); _stencil_val(p,0,-1,0);_stencil_val_a(gf.y,0,0,0);   }}end__stencil_is_face_y()
#line 403
_stencil_is_face_z(){
    {;; _stencil_val(alpha.z,0,0,0);_stencil_val(p,0,0,0); _stencil_val(p,0,0,-1);_stencil_val_a(gf.z,0,0,0);   }}end__stencil_is_face_z()}end_foreach_face_stencil();
  {
#line 403
foreach_face_generic(){is_face_x(){
    val(gf.x,0,0,0) = _const_fm.x*_const_a.x - val(alpha.x,0,0,0)*(val(p,0,0,0) - val(p,-1,0,0))/Delta;}end_is_face_x()
#line 403
is_face_y(){
    val(gf.y,0,0,0) = _const_fm.y*_const_a.y - val(alpha.y,0,0,0)*(val(p,0,0,0) - val(p,0,-1,0))/Delta;}end_is_face_y()
#line 403
is_face_z(){
    val(gf.z,0,0,0) = _const_fm.z*_const_a.z - val(alpha.z,0,0,0)*(val(p,0,0,0) - val(p,0,0,-1))/Delta;}end_is_face_z()}end_foreach_face_generic();}}else if(!is_constant(fm.x) && !is_constant(a.x) && is_constant(alpha.x)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX],_constant[alpha.z.i-_NVARMAX]};NOT_UNUSED(_const_alpha);
  
#line 403
foreach_face_stencil(1,{(NonLocal[]){{"p","scalar",(void *)&p,NULL,0},{"_const_alpha","_coord",(void *)&_const_alpha,NULL,0},{"a","vector",(void *)&a,NULL,0},{"fm","vector",(void *)&fm,NULL,0},{"gf","vector",(void *)&gf,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n// #line 403 \"/home/lisergey/basilisk/src/navier-stokes/centered.h\"\n{is_face_x(){\n    val_out_(gf.x,0,0,0) = val(fm.x,0,0,0)*val(a.x,0,0,0) - _const_alpha.x*(val(p,0,0,0) - val(p,-1,0,0))/Delta;}end_is_face_x()\n// #line 403\nis_face_y(){\n    val_out_(gf.y,0,0,0) = val(fm.y,0,0,0)*val(a.y,0,0,0) - _const_alpha.y*(val(p,0,0,0) - val(p,0,-1,0))/Delta;}end_is_face_y()\n// #line 403\nis_face_z(){\n    val_out_(gf.z,0,0,0) = val(fm.z,0,0,0)*val(a.z,0,0,0) - _const_alpha.z*(val(p,0,0,0) - val(p,0,0,-1))/Delta;}end_is_face_z()}"}){_stencil_is_face_x(){
    { _stencil_val(fm.x,0,0,0);_stencil_val(a.x,0,0,0);;_stencil_val(p,0,0,0); _stencil_val(p,-1,0,0);_stencil_val_a(gf.x,0,0,0);   }}end__stencil_is_face_x()
#line 403
_stencil_is_face_y(){
    { _stencil_val(fm.y,0,0,0);_stencil_val(a.y,0,0,0);;_stencil_val(p,0,0,0); _stencil_val(p,0,-1,0);_stencil_val_a(gf.y,0,0,0);   }}end__stencil_is_face_y()
#line 403
_stencil_is_face_z(){
    { _stencil_val(fm.z,0,0,0);_stencil_val(a.z,0,0,0);;_stencil_val(p,0,0,0); _stencil_val(p,0,0,-1);_stencil_val_a(gf.z,0,0,0);   }}end__stencil_is_face_z()}end_foreach_face_stencil();
  {
#line 403
foreach_face_generic(){is_face_x(){
    val(gf.x,0,0,0) = val(fm.x,0,0,0)*val(a.x,0,0,0) - _const_alpha.x*(val(p,0,0,0) - val(p,-1,0,0))/Delta;}end_is_face_x()
#line 403
is_face_y(){
    val(gf.y,0,0,0) = val(fm.y,0,0,0)*val(a.y,0,0,0) - _const_alpha.y*(val(p,0,0,0) - val(p,0,-1,0))/Delta;}end_is_face_y()
#line 403
is_face_z(){
    val(gf.z,0,0,0) = val(fm.z,0,0,0)*val(a.z,0,0,0) - _const_alpha.z*(val(p,0,0,0) - val(p,0,0,-1))/Delta;}end_is_face_z()}end_foreach_face_generic();}}else if(is_constant(fm.x) && !is_constant(a.x) && is_constant(alpha.x)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX],_constant[alpha.z.i-_NVARMAX]};NOT_UNUSED(_const_alpha);
  
#line 403
foreach_face_stencil(1,{(NonLocal[]){{"p","scalar",(void *)&p,NULL,0},{"_const_alpha","_coord",(void *)&_const_alpha,NULL,0},{"a","vector",(void *)&a,NULL,0},{"_const_fm","_coord",(void *)&_const_fm,NULL,0},{"gf","vector",(void *)&gf,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n// #line 403 \"/home/lisergey/basilisk/src/navier-stokes/centered.h\"\n{is_face_x(){\n    val_out_(gf.x,0,0,0) = _const_fm.x*val(a.x,0,0,0) - _const_alpha.x*(val(p,0,0,0) - val(p,-1,0,0))/Delta;}end_is_face_x()\n// #line 403\nis_face_y(){\n    val_out_(gf.y,0,0,0) = _const_fm.y*val(a.y,0,0,0) - _const_alpha.y*(val(p,0,0,0) - val(p,0,-1,0))/Delta;}end_is_face_y()\n// #line 403\nis_face_z(){\n    val_out_(gf.z,0,0,0) = _const_fm.z*val(a.z,0,0,0) - _const_alpha.z*(val(p,0,0,0) - val(p,0,0,-1))/Delta;}end_is_face_z()}"}){_stencil_is_face_x(){
    {;_stencil_val(a.x,0,0,0);;_stencil_val(p,0,0,0); _stencil_val(p,-1,0,0);_stencil_val_a(gf.x,0,0,0);   }}end__stencil_is_face_x()
#line 403
_stencil_is_face_y(){
    {;_stencil_val(a.y,0,0,0);;_stencil_val(p,0,0,0); _stencil_val(p,0,-1,0);_stencil_val_a(gf.y,0,0,0);   }}end__stencil_is_face_y()
#line 403
_stencil_is_face_z(){
    {;_stencil_val(a.z,0,0,0);;_stencil_val(p,0,0,0); _stencil_val(p,0,0,-1);_stencil_val_a(gf.z,0,0,0);   }}end__stencil_is_face_z()}end_foreach_face_stencil();
  {
#line 403
foreach_face_generic(){is_face_x(){
    val(gf.x,0,0,0) = _const_fm.x*val(a.x,0,0,0) - _const_alpha.x*(val(p,0,0,0) - val(p,-1,0,0))/Delta;}end_is_face_x()
#line 403
is_face_y(){
    val(gf.y,0,0,0) = _const_fm.y*val(a.y,0,0,0) - _const_alpha.y*(val(p,0,0,0) - val(p,0,-1,0))/Delta;}end_is_face_y()
#line 403
is_face_z(){
    val(gf.z,0,0,0) = _const_fm.z*val(a.z,0,0,0) - _const_alpha.z*(val(p,0,0,0) - val(p,0,0,-1))/Delta;}end_is_face_z()}end_foreach_face_generic();}}else if(!is_constant(fm.x) && is_constant(a.x) && is_constant(alpha.x)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX],_constant[a.z.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX],_constant[alpha.z.i-_NVARMAX]};NOT_UNUSED(_const_alpha);
  
#line 403
foreach_face_stencil(1,{(NonLocal[]){{"p","scalar",(void *)&p,NULL,0},{"_const_alpha","_coord",(void *)&_const_alpha,NULL,0},{"_const_a","_coord",(void *)&_const_a,NULL,0},{"fm","vector",(void *)&fm,NULL,0},{"gf","vector",(void *)&gf,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n// #line 403 \"/home/lisergey/basilisk/src/navier-stokes/centered.h\"\n{is_face_x(){\n    val_out_(gf.x,0,0,0) = val(fm.x,0,0,0)*_const_a.x - _const_alpha.x*(val(p,0,0,0) - val(p,-1,0,0))/Delta;}end_is_face_x()\n// #line 403\nis_face_y(){\n    val_out_(gf.y,0,0,0) = val(fm.y,0,0,0)*_const_a.y - _const_alpha.y*(val(p,0,0,0) - val(p,0,-1,0))/Delta;}end_is_face_y()\n// #line 403\nis_face_z(){\n    val_out_(gf.z,0,0,0) = val(fm.z,0,0,0)*_const_a.z - _const_alpha.z*(val(p,0,0,0) - val(p,0,0,-1))/Delta;}end_is_face_z()}"}){_stencil_is_face_x(){
    { _stencil_val(fm.x,0,0,0);;;_stencil_val(p,0,0,0); _stencil_val(p,-1,0,0);_stencil_val_a(gf.x,0,0,0);   }}end__stencil_is_face_x()
#line 403
_stencil_is_face_y(){
    { _stencil_val(fm.y,0,0,0);;;_stencil_val(p,0,0,0); _stencil_val(p,0,-1,0);_stencil_val_a(gf.y,0,0,0);   }}end__stencil_is_face_y()
#line 403
_stencil_is_face_z(){
    { _stencil_val(fm.z,0,0,0);;;_stencil_val(p,0,0,0); _stencil_val(p,0,0,-1);_stencil_val_a(gf.z,0,0,0);   }}end__stencil_is_face_z()}end_foreach_face_stencil();
  {
#line 403
foreach_face_generic(){is_face_x(){
    val(gf.x,0,0,0) = val(fm.x,0,0,0)*_const_a.x - _const_alpha.x*(val(p,0,0,0) - val(p,-1,0,0))/Delta;}end_is_face_x()
#line 403
is_face_y(){
    val(gf.y,0,0,0) = val(fm.y,0,0,0)*_const_a.y - _const_alpha.y*(val(p,0,0,0) - val(p,0,-1,0))/Delta;}end_is_face_y()
#line 403
is_face_z(){
    val(gf.z,0,0,0) = val(fm.z,0,0,0)*_const_a.z - _const_alpha.z*(val(p,0,0,0) - val(p,0,0,-1))/Delta;}end_is_face_z()}end_foreach_face_generic();}}else {_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX],_constant[a.z.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX],_constant[alpha.z.i-_NVARMAX]};NOT_UNUSED(_const_alpha);
  
#line 403
foreach_face_stencil(1,{(NonLocal[]){{"p","scalar",(void *)&p,NULL,0},{"_const_alpha","_coord",(void *)&_const_alpha,NULL,0},{"_const_a","_coord",(void *)&_const_a,NULL,0},{"_const_fm","_coord",(void *)&_const_fm,NULL,0},{"gf","vector",(void *)&gf,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n// #line 403 \"/home/lisergey/basilisk/src/navier-stokes/centered.h\"\n{is_face_x(){\n    val_out_(gf.x,0,0,0) = _const_fm.x*_const_a.x - _const_alpha.x*(val(p,0,0,0) - val(p,-1,0,0))/Delta;}end_is_face_x()\n// #line 403\nis_face_y(){\n    val_out_(gf.y,0,0,0) = _const_fm.y*_const_a.y - _const_alpha.y*(val(p,0,0,0) - val(p,0,-1,0))/Delta;}end_is_face_y()\n// #line 403\nis_face_z(){\n    val_out_(gf.z,0,0,0) = _const_fm.z*_const_a.z - _const_alpha.z*(val(p,0,0,0) - val(p,0,0,-1))/Delta;}end_is_face_z()}"}){_stencil_is_face_x(){
    {;;;_stencil_val(p,0,0,0); _stencil_val(p,-1,0,0);_stencil_val_a(gf.x,0,0,0);   }}end__stencil_is_face_x()
#line 403
_stencil_is_face_y(){
    {;;;_stencil_val(p,0,0,0); _stencil_val(p,0,-1,0);_stencil_val_a(gf.y,0,0,0);   }}end__stencil_is_face_y()
#line 403
_stencil_is_face_z(){
    {;;;_stencil_val(p,0,0,0); _stencil_val(p,0,0,-1);_stencil_val_a(gf.z,0,0,0);   }}end__stencil_is_face_z()}end_foreach_face_stencil();
  {
#line 403
foreach_face_generic(){is_face_x(){
    val(gf.x,0,0,0) = _const_fm.x*_const_a.x - _const_alpha.x*(val(p,0,0,0) - val(p,-1,0,0))/Delta;}end_is_face_x()
#line 403
is_face_y(){
    val(gf.y,0,0,0) = _const_fm.y*_const_a.y - _const_alpha.y*(val(p,0,0,0) - val(p,0,-1,0))/Delta;}end_is_face_y()
#line 403
is_face_z(){
    val(gf.z,0,0,0) = _const_fm.z*_const_a.z - _const_alpha.z*(val(p,0,0,0) - val(p,0,0,-1))/Delta;}end_is_face_z()}end_foreach_face_generic();}}





  trash (((vector[]){g,{{-1},{-1},{-1}}}));
  if(!is_constant(fm.x)){
  
#line 411
foreach_stencil(1,{(NonLocal[]){{"fm","vector",(void *)&fm,NULL,0},{"gf","vector",(void *)&gf,NULL,0},{"g","vector",(void *)&g,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n    \n// #line 412 \"/home/lisergey/basilisk/src/navier-stokes/centered.h\"\n{\n      val_out_(g.x,0,0,0) = (val(gf.x,0,0,0) + val(gf.x,1,0,0))/(val(fm.x,0,0,0) + val(fm.x,1,0,0) + 0.);\n      \n// #line 413\nval_out_(g.y,0,0,0) = (val(gf.y,0,0,0) + val(gf.y,0,1,0))/(val(fm.y,0,0,0) + val(fm.y,0,1,0) + 0.);\n      \n// #line 413\nval_out_(g.z,0,0,0) = (val(gf.z,0,0,0) + val(gf.z,0,0,1))/(val(fm.z,0,0,0) + val(fm.z,0,0,1) + 0.);}"})
    {
      {_stencil_val(gf.x,0,0,0); _stencil_val(gf.x,1,0,0);_stencil_val(fm.x,0,0,0); _stencil_val(fm.x,1,0,0);_stencil_val_a(g.x,0,0,0);      }
      
#line 413
{_stencil_val(gf.y,0,0,0); _stencil_val(gf.y,0,1,0);_stencil_val(fm.y,0,0,0); _stencil_val(fm.y,0,1,0);_stencil_val_a(g.y,0,0,0);      }
      
#line 413
{_stencil_val(gf.z,0,0,0); _stencil_val(gf.z,0,0,1);_stencil_val(fm.z,0,0,0); _stencil_val(fm.z,0,0,1);_stencil_val_a(g.z,0,0,0);      }}end_foreach_stencil();{
#line 411
foreach()
    {
      val(g.x,0,0,0) = (val(gf.x,0,0,0) + val(gf.x,1,0,0))/(val(fm.x,0,0,0) + val(fm.x,1,0,0) + 0.);
      
#line 413
val(g.y,0,0,0) = (val(gf.y,0,0,0) + val(gf.y,0,1,0))/(val(fm.y,0,0,0) + val(fm.y,0,1,0) + 0.);
      
#line 413
val(g.z,0,0,0) = (val(gf.z,0,0,0) + val(gf.z,0,0,1))/(val(fm.z,0,0,0) + val(fm.z,0,0,1) + 0.);}end_foreach();}}else {_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);
  
#line 411
foreach_stencil(1,{(NonLocal[]){{"_const_fm","_coord",(void *)&_const_fm,NULL,0},{"gf","vector",(void *)&gf,NULL,0},{"g","vector",(void *)&g,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n    \n// #line 412 \"/home/lisergey/basilisk/src/navier-stokes/centered.h\"\n{\n      val_out_(g.x,0,0,0) = (val(gf.x,0,0,0) + val(gf.x,1,0,0))/(_const_fm.x + _const_fm.x + 0.);\n      \n// #line 413\nval_out_(g.y,0,0,0) = (val(gf.y,0,0,0) + val(gf.y,0,1,0))/(_const_fm.y + _const_fm.y + 0.);\n      \n// #line 413\nval_out_(g.z,0,0,0) = (val(gf.z,0,0,0) + val(gf.z,0,0,1))/(_const_fm.z + _const_fm.z + 0.);}"})
    {
      {_stencil_val(gf.x,0,0,0); _stencil_val(gf.x,1,0,0);;;_stencil_val_a(g.x,0,0,0);      }
      
#line 413
{_stencil_val(gf.y,0,0,0); _stencil_val(gf.y,0,1,0);;;_stencil_val_a(g.y,0,0,0);      }
      
#line 413
{_stencil_val(gf.z,0,0,0); _stencil_val(gf.z,0,0,1);;;_stencil_val_a(g.z,0,0,0);      }}end_foreach_stencil();
  {
#line 411
foreach()
    {
      val(g.x,0,0,0) = (val(gf.x,0,0,0) + val(gf.x,1,0,0))/(_const_fm.x + _const_fm.x + 0.);
      
#line 413
val(g.y,0,0,0) = (val(gf.y,0,0,0) + val(gf.y,0,1,0))/(_const_fm.y + _const_fm.y + 0.);
      
#line 413
val(g.z,0,0,0) = (val(gf.z,0,0,0) + val(gf.z,0,0,1))/(_const_fm.z + _const_fm.z + 0.);}end_foreach();}}delete((scalar*)((vector[]){gf,{{-1},{-1},{-1}}}));
}






static int projection_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i++)!=0;*ip=i;*tp=t;return ret;}







#line 421
      static int projection(const int i,const double t,Event *_ev){tracing("projection","/home/lisergey/basilisk/src/navier-stokes/centered.h",421);
{
  mgp = project (uf, p, alpha, dt, mgp.nrelax);
  centered_gradient (p, g);




  correction (dt);
}{end_tracing("projection","/home/lisergey/basilisk/src/navier-stokes/centered.h",430);return 0;}end_tracing("projection","/home/lisergey/basilisk/src/navier-stokes/centered.h",430);}





static int end_timestep_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i++)!=0;*ip=i;*tp=t;return ret;}






#line 436
static int end_timestep(const int i,const double t,Event *_ev){;return 0;}









static int adapt_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i++)!=0;*ip=i;*tp=t;return ret;}
#line 446 "/home/lisergey/basilisk/src/navier-stokes/centered.h"
      static int adapt(const int i,const double t,Event *_ev){tracing("adapt","/home/lisergey/basilisk/src/navier-stokes/centered.h",446); {






  event ("properties");
}{end_tracing("adapt","/home/lisergey/basilisk/src/navier-stokes/centered.h",454);return 0;}end_tracing("adapt","/home/lisergey/basilisk/src/navier-stokes/centered.h",454);}
#line 10 "cylinder.c"
#line 1 "lambda2.h"
#line 1 "/home/lisergey/basilisk/src/lambda2.h"
static void eigsrt (double d[3],
      double v[3][3])
{
  int k, j, i;
  double p;

  for (i = 0; i < 3 - 1; i++) {
    p = d[k = i];

    for (j = i + 1; j < 3; j++)
      if (d[j] >= p)
 p = d[k = j];
    if (k != i) {
      d[k] = d[i];
      d[i] = p;
      for (j = 0; j < 3; j++) {
 p = v[j][i];
 v[j][i] = v[j][k];
 v[j][k] = p;
      }
    }
  }
}
#line 37 "/home/lisergey/basilisk/src/lambda2.h"
void eigenvalues (double a[3][3],
    double d[3],
    double v[3][3])
{
  int j, iq, ip, i;
  double tresh, theta, tau, t, sm, s, h, g, c, b[3], z[3];

  for (ip = 0; ip < 3; ip++) {
    for (iq = 0; iq < 3; iq++)
      v[ip][iq] = 0.0;
    v[ip][ip] = 1.0;
  }

  for (ip = 0; ip < 3; ip++) {
    b[ip] = d[ip] = a[ip][ip];
    z[ip] = 0.0;
  }

  for (i = 1; i <= 50; i++) {
    sm = 0.0;
    for (ip = 0; ip < 3 - 1; ip++) {
      for (iq = ip + 1; iq < 3; iq++)
 sm += fabs (a[ip][iq]);
    }
    if (sm == 0.0) {
      eigsrt (d, v);
      return;
    }
    if (i < 4)
      tresh = 0.2*sm/(3*3);
    else
      tresh = 0.0;
    for (ip = 0; ip < 3 - 1; ip++) {
      for (iq = ip + 1; iq < 3; iq++) {
 g = 100.0*fabs (a[ip][iq]);
 if (i > 4 && fabs(d[ip]) + g == fabs(d[ip]) &&
     fabs(d[iq]) + g == fabs(d[iq]))
   a[ip][iq] = 0.0;
 else if (fabs (a[ip][iq]) > tresh) {
   h = d[iq] - d[ip];
   if (fabs(h) + g == fabs(h))
     t = a[ip][iq]/h;
   else {
     theta = 0.5*h/a[ip][iq];
     t = 1.0/(fabs (theta) + sqrt (1.0 + theta*theta));
     if (theta < 0.0) t = -t;
   }
   c = 1.0/sqrt (1 + t*t);
   s = t*c;
   tau = s/(1.0 + c);
   h = t*a[ip][iq];
   z[ip] -= h;
   z[iq] += h;
   d[ip] -= h;
   d[iq] += h;
   a[ip][iq] = 0.0;
   for (j = 0; j <= ip - 1; j++)
     { g=a[j][ip];h=a[j][iq];a[j][ip]=g-s*(h+g*tau);a[j][iq]=h+s*(g-h*tau);};
   for (j = ip + 1; j <= iq - 1; j++)
     { g=a[ip][j];h=a[j][iq];a[ip][j]=g-s*(h+g*tau);a[j][iq]=h+s*(g-h*tau);};
   for (j = iq + 1; j < 3; j++)
     { g=a[ip][j];h=a[iq][j];a[ip][j]=g-s*(h+g*tau);a[iq][j]=h+s*(g-h*tau);};
   for (j = 0; j < 3; j++)
     { g=v[j][ip];h=v[j][iq];v[j][ip]=g-s*(h+g*tau);v[j][iq]=h+s*(g-h*tau);};
 }
      }
    }
    for (ip = 0; ip < 3; ip++) {
      b[ip] += z[ip];
      d[ip] = b[ip];
      z[ip] = 0.0;
    }
  }

  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++)
      fprintf (ferr, "%10.3g ", a[i][j]);
    fprintf (ferr, "\n");
  }
  if (!(false)) qassert ("/home/lisergey/basilisk/src/lambda2.h", 116, "false");
}

void lambda2 (const vector u, scalar l2)
{
  foreach_stencil(1,{(NonLocal[]){{"l2","scalar",(void *)&l2,NULL,0},{"ferr","not implemented yet",(void *)ferr,NULL,1},{"u","vector",(void *)&u,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"\n@ def not_mpi_compatible()\ndo {\n  if (npe() > 1) {\n    fprintf (ferr, \"%s() is not compatible with MPI (yet)\\n\", __func__);\n    exit (1);\n  }\n} while(0)\n@\n@ define system(command) (pid() == 0 ? system(command) : 0)\n@else\n@ define qstderr() stderr\n@ define qstdout() stdout\n@ define ferr stderr\n@ define fout stdout\n@ define not_mpi_compatible()\n@endif\n\n\n\n\n// #line 93 \"/home/lisergey/basilisk/src/common.h\"\nstatic inline void qassert (const char * file, int line, const char * cond) {\n  fprintf (ferr, \"%s:%d: Assertion `%s' failed.\\n\", file, line, cond);\n  abort();\n}\n// # 10 \"cylinder.c\" 2\n// # 1 \"lambda2.h\" 1\n// # 1 \"/home/lisergey/basilisk/src/lambda2.h\"\nstatic void eigsrt (real d[3],\n      real v[3][3])\n{\n  int k, j, i;\n  real p;\n\n  for (i = 0; i < 3 - 1; i++) {\n    p = d[k = i];\n\n    for (j = i + 1; j < 3; j++)\n      if (d[j] >= p)\n p = d[k = j];\n    if (k != i) {\n      d[k] = d[i];\n      d[i] = p;\n      for (j = 0; j < 3; j++) {\n p = v[j][i];\n v[j][i] = v[j][k];\n v[j][k] = p;\n      }\n    }\n  }\n}\n// # 37 \"/home/lisergey/basilisk/src/lambda2.h\"\nvoid eigenvalues (real a[3][3],\n    real d[3],\n    real v[3][3])\n{\n  int j, iq, ip, i;\n  real tresh, theta, tau, t, sm, s, h, g, c, b[3], z[3];\n\n  for (ip = 0; ip < 3; ip++) {\n    for (iq = 0; iq < 3; iq++)\n      v[ip][iq] = 0.0;\n    v[ip][ip] = 1.0;\n  }\n\n  for (ip = 0; ip < 3; ip++) {\n    b[ip] = d[ip] = a[ip][ip];\n    z[ip] = 0.0;\n  }\n\n  for (i = 1; i <= 50; i++) {\n    sm = 0.0;\n    for (ip = 0; ip < 3 - 1; ip++) {\n      for (iq = ip + 1; iq < 3; iq++)\n sm += fabs (a[ip][iq]);\n    }\n    if (sm == 0.0) {\n      eigsrt (d, v);\n      return;\n    }\n    if (i < 4)\n      tresh = 0.2*sm/(3*3);\n    else\n      tresh = 0.0;\n    for (ip = 0; ip < 3 - 1; ip++) {\n      for (iq = ip + 1; iq < 3; iq++) {\n g = 100.0*fabs (a[ip][iq]);\n if (i > 4 && fabs(d[ip]) + g == fabs(d[ip]) &&\n     fabs(d[iq]) + g == fabs(d[iq]))\n   a[ip][iq] = 0.0;\n else if (fabs (a[ip][iq]) > tresh) {\n   h = d[iq] - d[ip];\n   if (fabs(h) + g == fabs(h))\n     t = a[ip][iq]/h;\n   else {\n     theta = 0.5*h/a[ip][iq];\n     t = 1.0/(fabs (theta) + sqrt (1.0 + theta*theta));\n     if (theta < 0.0) t = -t;\n   }\n   c = 1.0/sqrt (1 + t*t);\n   s = t*c;\n   tau = s/(1.0 + c);\n   h = t*a[ip][iq];\n   z[ip] -= h;\n   z[iq] += h;\n   d[ip] -= h;\n   d[iq] += h;\n   a[ip][iq] = 0.0;\n   for (j = 0; j <= ip - 1; j++)\n     { g=a[j][ip];h=a[j][iq];a[j][ip]=g-s*(h+g*tau);a[j][iq]=h+s*(g-h*tau);};\n   for (j = ip + 1; j <= iq - 1; j++)\n     { g=a[ip][j];h=a[j][iq];a[ip][j]=g-s*(h+g*tau);a[j][iq]=h+s*(g-h*tau);};\n   for (j = iq + 1; j < 3; j++)\n     { g=a[ip][j];h=a[iq][j];a[ip][j]=g-s*(h+g*tau);a[iq][j]=h+s*(g-h*tau);};\n   for (j = 0; j < 3; j++)\n     { g=v[j][ip];h=v[j][iq];v[j][ip]=g-s*(h+g*tau);v[j][iq]=h+s*(g-h*tau);};\n }\n      }\n    }\n    for (ip = 0; ip < 3; ip++) {\n      b[ip] += z[ip];\n      d[ip] = b[ip];\n      z[ip] = 0.0;\n    }\n  }\n\n  for (i = 0; i < 3; i++) {\n    for (j = 0; j < 3; j++)\n      fprintf (ferr, \"%10.3g \", a[i][j]);\n    fprintf (ferr, \"\\n\");\n  }\n  if (!(false)) qassert (\"/home/lisergey/basilisk/src/lambda2.h\", 116, \"false\");\n}"," \n// #line 121 \"/home/lisergey/basilisk/src/lambda2.h\"\n{\n    real JJ[3][3];\n    scalar s = u.x;\n    int i = 0;\n    \n      JJ[0][i++] = ((val(s,1,0,0) - val(s,-1,0,0))/(2.*Delta));\n      \n// #line 126\nJJ[0][i++] = ((val(s,0,1,0) - val(s,0,-1,0))/(2.*Delta));\n      \n// #line 126\nJJ[0][i++] = ((val(s,0,0,1) - val(s,0,0,-1))/(2.*Delta));\n    s = u.y; i = 0;\n    \n      JJ[1][i++] = ((val(s,1,0,0) - val(s,-1,0,0))/(2.*Delta));\n      \n// #line 129\nJJ[1][i++] = ((val(s,0,1,0) - val(s,0,-1,0))/(2.*Delta));\n      \n// #line 129\nJJ[1][i++] = ((val(s,0,0,1) - val(s,0,0,-1))/(2.*Delta));\n    s = u.z; i = 0;\n    \n      JJ[2][i++] = ((val(s,1,0,0) - val(s,-1,0,0))/(2.*Delta));\n      \n// #line 132\nJJ[2][i++] = ((val(s,0,1,0) - val(s,0,-1,0))/(2.*Delta));\n      \n// #line 132\nJJ[2][i++] = ((val(s,0,0,1) - val(s,0,0,-1))/(2.*Delta));\n    real S2O2[3][3];\n    for (int i = 0; i < 3; i++)\n      for (int j = 0; j < 3; j++) {\n S2O2[i][j] = 0.;\n for (int k = 0; k < 3; k++)\n   S2O2[i][j] += JJ[i][k]*JJ[k][j] + JJ[k][i]*JJ[j][k];\n      }\n    real lambda[3], ev[3][3];\n    eigenvalues (S2O2, lambda, ev);\n    val_out_(l2,0,0,0) = lambda[1]/2.;\n  }"}) { 
    
    scalar s = u.x;   
    
    
      {_stencil_val(s,1,0,0); _stencil_val(s,-1,0,0);   }
      
#line 126
{_stencil_val(s,0,1,0); _stencil_val(s,0,-1,0);   }
      
#line 126
{_stencil_val(s,0,0,1); _stencil_val(s,0,0,-1);   }
    s = u.y;   
    
      {_stencil_val(s,1,0,0); _stencil_val(s,-1,0,0);   }
      
#line 129
{_stencil_val(s,0,1,0); _stencil_val(s,0,-1,0);   }
      
#line 129
{_stencil_val(s,0,0,1); _stencil_val(s,0,0,-1);   }
    s = u.z;   
    
      {_stencil_val(s,1,0,0); _stencil_val(s,-1,0,0);   }
      
#line 132
{_stencil_val(s,0,1,0); _stencil_val(s,0,-1,0);   }
      
#line 132
{_stencil_val(s,0,0,1); _stencil_val(s,0,0,-1);   }   
            
            
       
       
     
      
       
       
    
    _stencil_val_a(l2,0,0,0);  
  }end_foreach_stencil();
  {
#line 121
foreach() {
    double JJ[3][3];
    scalar s = u.x;
    int i = 0;
    
      JJ[0][i++] = ((val(s,1,0,0) - val(s,-1,0,0))/(2.*Delta));
      
#line 126
JJ[0][i++] = ((val(s,0,1,0) - val(s,0,-1,0))/(2.*Delta));
      
#line 126
JJ[0][i++] = ((val(s,0,0,1) - val(s,0,0,-1))/(2.*Delta));
    s = u.y; i = 0;
    
      JJ[1][i++] = ((val(s,1,0,0) - val(s,-1,0,0))/(2.*Delta));
      
#line 129
JJ[1][i++] = ((val(s,0,1,0) - val(s,0,-1,0))/(2.*Delta));
      
#line 129
JJ[1][i++] = ((val(s,0,0,1) - val(s,0,0,-1))/(2.*Delta));
    s = u.z; i = 0;
    
      JJ[2][i++] = ((val(s,1,0,0) - val(s,-1,0,0))/(2.*Delta));
      
#line 132
JJ[2][i++] = ((val(s,0,1,0) - val(s,0,-1,0))/(2.*Delta));
      
#line 132
JJ[2][i++] = ((val(s,0,0,1) - val(s,0,0,-1))/(2.*Delta));
    double S2O2[3][3];
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++) {
 S2O2[i][j] = 0.;
 for (int k = 0; k < 3; k++)
   S2O2[i][j] += JJ[i][k]*JJ[k][j] + JJ[k][i]*JJ[j][k];
      }
    double lambda[3], ev[3][3];
    eigenvalues (S2O2, lambda, ev);
    val(l2,0,0,0) = lambda[1]/2.;
  }end_foreach();}
}
#line 11 "cylinder.c"
#line 1 "output_xdmf.h"
#line 1 "./output_xdmf.h"
static int output_xdmf(double t, scalar *list, vector *vlist,
                       int cond(double, double, double, double),
                       const char *path) {
  float *xyz, *attr;
  long j, nattr, nvect, ncell, ncell_total, nsize, offset;
  char xyz_path[FILENAME_MAX], attr_path[FILENAME_MAX], xdmf_path[FILENAME_MAX],
      *vname, *xyz_base, *attr_base;
  FILE *file;
  MPI_File mpi_file;
  const int shift[8][3] = {
      {0, 0, 0}, {0, 0, 1}, {0, 1, 1}, {0, 1, 0},
      {1, 0, 0}, {1, 0, 1}, {1, 1, 1}, {1, 1, 0},
  };

  snprintf(xyz_path, sizeof xyz_path, "%s.xyz.raw", path);
  snprintf(attr_path, sizeof attr_path, "%s.attr.raw", path);
  snprintf(xdmf_path, sizeof xdmf_path, "%s.xdmf2", path);

  xyz_base = xyz_path;
  attr_base = attr_path;
  for (j = 0; xyz_path[j] != '\0'; j++) {
    if (xyz_path[j] == '/' && xyz_path[j + 1] != '\0') {
      xyz_base = &xyz_path[j + 1];
      attr_base = &attr_path[j + 1];
    }
  }

  nsize = 0;
  ncell = 0;
  j = 0;
  xyz = NULL;
  {foreach_cell() if (is_local(cell) && is_leaf(cell) &&
                     (!cond || cond(x, y, z, Delta))) {
    long i;
    ncell++;
    if (ncell >= nsize) {
      nsize = 2 * nsize + 1;
      if ((xyz = prealloc(xyz, 8 * 3 * nsize * sizeof *xyz,__func__,__FILE__,__LINE__)) == NULL) {
        fprintf(ferr, "%s:%d: realloc failed\n", "./output_xdmf.h", 39);
        return 1;
      }
    }
    for (i = 0; i < 8; i++) {
      xyz[j++] = x + Delta * (shift[i][0] - 0.5);
      xyz[j++] = y + Delta * (shift[i][1] - 0.5);
      xyz[j++] = z + Delta * (shift[i][2] - 0.5);
    }
  }end_foreach_cell();}

  if ((file = fopen(xyz_path, "w")) == NULL) {
    fprintf(ferr, "%s:%d: fail to open '%s'\n", "./output_xdmf.h", 51, xyz_path);
    return 1;
  }

  MPI_Exscan(&ncell, &offset, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
  if (pid() == 0)
    offset = 0;
  MPI_File_open(MPI_COMM_WORLD, xyz_path, MPI_MODE_CREATE | MPI_MODE_WRONLY,
                MPI_INFO_NULL, &mpi_file);
  MPI_File_write_at_all(mpi_file, 3 * 8 * offset * sizeof *xyz, xyz,
                        3 * 8 * ncell * sizeof *xyz, MPI_BYTE,
                        MPI_STATUS_IGNORE);
  pfree(xyz,__func__,__FILE__,__LINE__);
  MPI_File_close(&mpi_file);

  nattr = list_len(list);
  nvect = vectors_len(vlist);
  if ((attr = pmalloc((nattr + 3 * nvect) * ncell * sizeof *attr,__func__,__FILE__,__LINE__)) == NULL) {
    fprintf(ferr, "%s:%d: malloc failed\n", "./output_xdmf.h", 69);
    return 1;
  }
  j = 0;
  {foreach_cell() if (is_local(cell) && is_leaf(cell) &&
                     (!cond || cond(x, y, z, Delta))) {
    {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      attr[j++] = val(s,0,0,0);}}
    {vector*_i=(vector*)( vlist);if(_i)for(vector v=*_i;(&v)->x.i>=0;v=*++_i){ {
      attr[j++] = val(v.x,0,0,0);
      attr[j++] = val(v.y,0,0,0);
      attr[j++] = val(v.z,0,0,0);
    }}}
  }end_foreach_cell();}
  if (!(j == (nattr + 3 * nvect) * ncell)) qassert ("./output_xdmf.h", 83, "j == (nattr + 3 * nvect) * ncell");
  MPI_File_open(MPI_COMM_WORLD, attr_path, MPI_MODE_CREATE | MPI_MODE_WRONLY,
                MPI_INFO_NULL, &mpi_file);
  MPI_File_write_at_all(mpi_file, (nattr + 3 * nvect) * offset * sizeof *attr,
                        attr, (nattr + 3 * nvect) * ncell * sizeof *attr,
                        MPI_BYTE, MPI_STATUS_IGNORE);
  pfree(attr,__func__,__FILE__,__LINE__);
  MPI_File_close(&mpi_file);

  if (pid() == npe() - 1) {
    ncell_total = offset + ncell;
    if ((file = fopen(xdmf_path, "w")) == NULL) {
      fprintf(ferr, "%s:%d: fail to open '%s'\n", "./output_xdmf.h", 95,
              xdmf_path);
      return 1;
    }
    fprintf(file,
            "<Xdmf\n"
            "    Version=\"2\">\n"
            "  <Domain>\n"
            "    <Grid>\n"
            "      <Time\n"
            "          Value=\"%.16e\"/>\n"
            "      <Topology\n"
            "          TopologyType=\"Hexahedron\"\n"
            "          Dimensions=\"%ld\"/>\n"
            "      <Geometry>\n"
            "        <DataItem\n"
            "            Dimensions=\"%ld 3\"\n"
            "            Format=\"Binary\">\n"
            "          %s\n"
            "        </DataItem>\n"
            "      </Geometry>\n",
            t, ncell_total, 8 * ncell_total, xyz_base);
    j = 0;
    {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      fprintf(file,
              "      <Attribute\n"
              "          Name=\"%s\"\n"
              "          Center=\"Cell\">\n"
              "        <DataItem\n"
              "            ItemType=\"HyperSlab\"\n"
              "            Dimensions=\"%ld\"\n"
              "            Type=\"HyperSlab\">\n"
              "          <DataItem Dimensions=\"3 1\">\n"
              "            %ld %ld %ld\n"
              "          </DataItem>\n"
              "          <DataItem\n"
              "              Dimensions=\"%ld\"\n"
              "              Format=\"Binary\">\n"
              "            %s\n"
              "          </DataItem>\n"
              "         </DataItem>\n"
              "      </Attribute>\n",
              _attribute[s.i].name, ncell_total, j++, nattr + 3 * nvect, ncell_total,
              (nattr + 3 * nvect) * ncell_total, attr_base);}}
    {vector*_i=(vector*)( vlist);if(_i)for(vector v=*_i;(&v)->x.i>=0;v=*++_i){ {
      vname = pstrdup(_attribute[v.x.i].name,__func__,__FILE__,__LINE__);
      *strrchr(vname, '.') = '\0';
      fprintf(file,
              "      <Attribute\n"
              "          Name=\"%s\"\n"
              "          AttributeType=\"Vector\"\n"
              "          Center=\"Cell\">\n"
              "        <DataItem\n"
              "            ItemType=\"HyperSlab\"\n"
              "            Dimensions=\"%ld 3\"\n"
              "            Type=\"HyperSlab\">\n"
              "          <DataItem Dimensions=\"3 2\">\n"
              "            0 %ld\n"
              "            1 1\n"
              "            %ld 3\n"
              "          </DataItem>\n"
              "          <DataItem\n"
              "              Dimensions=\"%ld %ld\"\n"
              "              Format=\"Binary\">\n"
              "            %s\n"
              "          </DataItem>\n"
              "         </DataItem>\n"
              "      </Attribute>\n",
              vname, ncell_total, j, ncell_total, ncell_total,
              nattr + 3 * nvect, attr_base);
      pfree(vname,__func__,__FILE__,__LINE__);
      j += 3;
    }}}
    fprintf(file, "    </Grid>\n"
                  "  </Domain>\n"
                  "</Xdmf>\n");
    if (fclose(file) != 0) {
      fprintf(ferr, "%s:%d: error: fail to close '%s'\n", "./output_xdmf.h", 172,
              xdmf_path);
      return 1;
    }
  }
  return 0;
}
#line 12 "cylinder.c"
#line 79 "cylinder.c"
coord Force = {0};
scalar  cs={11};
vector  fs={{12},{13},{14}};
      static void embed_force3(scalar p, vector u, vector mu, coord *Fp,
                               coord *Fmu) {tracing("embed_force3","cylinder.c",82);
   {
    Fmu->x = Force.x;
    Fp->x = 0;
  } 
#line 84
{
    Fmu->y = Force.y;
    Fp->y = 0;
  } 
#line 84
{
    Fmu->z = Force.z;
    Fp->z = 0;
  }
end_tracing("embed_force3","cylinder.c",88);}
struct Cleanup {
  scalar c;
  vector s;
  double smin;
  bool opposite;
};
      static int fractions_cleanup(struct Cleanup u) {tracing("fractions_cleanup","cylinder.c",95); {end_tracing("fractions_cleanup","cylinder.c",95);return 0;} end_tracing("fractions_cleanup","cylinder.c",95);}
static int velocity_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i++)!=0;*ip=i;*tp=t;return ret;}

#line 96
      static int velocity(const int i,const double t,Event *_ev){tracing("velocity","cylinder.c",96); {
   Force.x = 0; Force.y = 0; Force.z = 0;
  foreach_stencil (1,{(NonLocal[]){{"Force","coord",(void *)&Force,NULL,0,'+'},{"cs","scalar",(void *)&cs,NULL,0},{"u","vector",(void *)&u,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},""," \n// #line 98 \"cylinder.c\"\n{\n     {\n      val_out_(u.x,0,0,0) = val(cs,0,0,0) * val(u.x,0,0,0);\n      Force.x += val(u.x,0,0,0) * (val(cs,0,0,0) - 1) * Delta * Delta * Delta;\n    } \n// #line 99\n{\n      val_out_(u.y,0,0,0) = val(cs,0,0,0) * val(u.y,0,0,0);\n      Force.y += val(u.y,0,0,0) * (val(cs,0,0,0) - 1) * Delta * Delta * Delta;\n    } \n// #line 99\n{\n      val_out_(u.z,0,0,0) = val(cs,0,0,0) * val(u.z,0,0,0);\n      Force.z += val(u.z,0,0,0) * (val(cs,0,0,0) - 1) * Delta * Delta * Delta;\n    }\n  }"}) {
     { 
_stencil_val(cs,0,0,0); _stencil_val(u.x,0,0,0);
      
#line 100
_stencil_val_a(u.x,0,0,0); 
_stencil_val(u.x,0,0,0);_stencil_val(cs,0,0,0);   
                
    
#line 102
} 
#line 99
{ 
_stencil_val(cs,0,0,0); _stencil_val(u.y,0,0,0);
      
#line 100
_stencil_val_a(u.y,0,0,0); 
_stencil_val(u.y,0,0,0);_stencil_val(cs,0,0,0);   
                
    
#line 102
} 
#line 99
{ 
_stencil_val(cs,0,0,0); _stencil_val(u.z,0,0,0);
      
#line 100
_stencil_val_a(u.z,0,0,0); 
_stencil_val(u.z,0,0,0);_stencil_val(cs,0,0,0);   
                
    
#line 102
}
  }end_foreach_stencil();
  
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel reduction(+ : Force)){
#line 98
foreach () {
     {
      val(u.x,0,0,0) = val(cs,0,0,0) * val(u.x,0,0,0);
      Force.x += val(u.x,0,0,0) * (val(cs,0,0,0) - 1) * Delta * Delta * Delta;
    } 
#line 99
{
      val(u.y,0,0,0) = val(cs,0,0,0) * val(u.y,0,0,0);
      Force.y += val(u.y,0,0,0) * (val(cs,0,0,0) - 1) * Delta * Delta * Delta;
    } 
#line 99
{
      val(u.z,0,0,0) = val(cs,0,0,0) * val(u.z,0,0,0);
      Force.z += val(u.z,0,0,0) * (val(cs,0,0,0) - 1) * Delta * Delta * Delta;
    }
  }end_foreach();mpi_all_reduce_array((double *)&Force,double,MPI_SUM,sizeof(Force)/sizeof(double));
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
   
#line 104
Force.x /= dt; Force.y /= dt; Force.z /= dt;
}{end_tracing("velocity","cylinder.c",105);return 0;}end_tracing("velocity","cylinder.c",105);}


static double dot3(const double *a, const double *b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

static void vorticity_vector(const vector u, vector omega) {
  if(!is_constant(cm) && !is_constant(fm.x)){
  
#line 113
foreach_stencil (1,{(NonLocal[]){{"omega","vector",(void *)&omega,NULL,0},{"u","vector",(void *)&u,NULL,0},{"fm","vector",(void *)&fm,NULL,0},{"cm","scalar",(void *)&cm,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"\n\n\n\n// #line 108 \"cylinder.c\"\nstatic real dot3(const real *a, const real *b) {\n  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];\n}"," \n// #line 113 \"cylinder.c\"\n{\n    real delta;\n    delta = (2. * val(cm,0,0,0) * Delta + 0.);\n    real fx[3] = {val(fm.x,1,0,0) - val(fm.x,0,0,0), val(fm.x,1,0,0), -val(fm.x,0,0,0)};\n    real fy[3] = {val(fm.y,0, 1,0) - val(fm.y,0,0,0), val(fm.y,0, 1,0), -val(fm.y,0,0,0)};\n    real fz[3] = {val(fm.z,0, 0, 1) - val(fm.z,0,0,0), val(fm.z,0, 0, 1), -val(fm.z,0,0,0)};\n    real xy[3] = {val(u.x,0,0,0), val(u.x,0, 1,0), val(u.x,0, -1,0)};\n    real xz[3] = {val(u.x,0,0,0), val(u.x,0, 0, 1), val(u.x,0, 0, -1)};\n    real yx[3] = {val(u.y,0,0,0), val(u.y,1,0,0), val(u.y,-1,0,0)};\n    real yz[3] = {val(u.y,0,0,0), val(u.y,0, 0, 1), val(u.y,0, 0, -1)};\n    real zx[3] = {val(u.z,0,0,0), val(u.z,1,0,0), val(u.z,-1,0,0)};\n    real zy[3] = {val(u.z,0,0,0), val(u.z,0, 1,0), val(u.z,0, -1,0)};\n    val_out_(omega.x,0,0,0) = (dot3(fy, zy) - dot3(fz, yz)) / delta;\n    val_out_(omega.y,0,0,0) = (dot3(fz, xz) - dot3(fx, zx)) / delta;\n    val_out_(omega.z,0,0,0) = (dot3(fx, yx) - dot3(fy, xy)) / delta;\n  }"}) { 
     
_stencil_val(cm,0,0,0);      
          
    
#line 116
_stencil_val(fm.x,0,0,0); _stencil_val(fm.x,1,0,0); _stencil_val(fm.x,0,0,0);_stencil_val(fm.x,1,0,0);     
    _stencil_val(fm.y,0,0,0); _stencil_val(fm.y,0, 1,0); _stencil_val(fm.y,0,0,0);_stencil_val(fm.y,0, 1,0);     
    _stencil_val(fm.z,0,0,0); _stencil_val(fm.z,0, 0, 1); _stencil_val(fm.z,0,0,0);_stencil_val(fm.z,0, 0, 1);   
     _stencil_val(u.x,0, -1,0); _stencil_val(u.x,0, 1,0);_stencil_val(u.x,0,0,0);   
     _stencil_val(u.x,0, 0, -1); _stencil_val(u.x,0, 0, 1);_stencil_val(u.x,0,0,0);   
     _stencil_val(u.y,-1,0,0); _stencil_val(u.y,1,0,0);_stencil_val(u.y,0,0,0);   
     _stencil_val(u.y,0, 0, -1); _stencil_val(u.y,0, 0, 1);_stencil_val(u.y,0,0,0);   
     _stencil_val(u.z,-1,0,0); _stencil_val(u.z,1,0,0);_stencil_val(u.z,0,0,0);   
     _stencil_val(u.z,0, -1,0); _stencil_val(u.z,0, 1,0);_stencil_val(u.z,0,0,0);
    _stencil_val_a(omega.x,0,0,0);        
    _stencil_val_a(omega.y,0,0,0);        
    _stencil_val_a(omega.z,0,0,0);        
  }end_foreach_stencil();{
#line 113
foreach () {
    double delta;
    delta = (2. * val(cm,0,0,0) * Delta + 0.);
    double fx[3] = {val(fm.x,1,0,0) - val(fm.x,0,0,0), val(fm.x,1,0,0), -val(fm.x,0,0,0)};
    double fy[3] = {val(fm.y,0, 1,0) - val(fm.y,0,0,0), val(fm.y,0, 1,0), -val(fm.y,0,0,0)};
    double fz[3] = {val(fm.z,0, 0, 1) - val(fm.z,0,0,0), val(fm.z,0, 0, 1), -val(fm.z,0,0,0)};
    double xy[3] = {val(u.x,0,0,0), val(u.x,0, 1,0), val(u.x,0, -1,0)};
    double xz[3] = {val(u.x,0,0,0), val(u.x,0, 0, 1), val(u.x,0, 0, -1)};
    double yx[3] = {val(u.y,0,0,0), val(u.y,1,0,0), val(u.y,-1,0,0)};
    double yz[3] = {val(u.y,0,0,0), val(u.y,0, 0, 1), val(u.y,0, 0, -1)};
    double zx[3] = {val(u.z,0,0,0), val(u.z,1,0,0), val(u.z,-1,0,0)};
    double zy[3] = {val(u.z,0,0,0), val(u.z,0, 1,0), val(u.z,0, -1,0)};
    val(omega.x,0,0,0) = (dot3(fy, zy) - dot3(fz, yz)) / delta;
    val(omega.y,0,0,0) = (dot3(fz, xz) - dot3(fx, zx)) / delta;
    val(omega.z,0,0,0) = (dot3(fx, yx) - dot3(fy, xy)) / delta;
  }end_foreach();}}else if(is_constant(cm) && !is_constant(fm.x)){double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);
  
#line 113
foreach_stencil (1,{(NonLocal[]){{"omega","vector",(void *)&omega,NULL,0},{"u","vector",(void *)&u,NULL,0},{"fm","vector",(void *)&fm,NULL,0},{"_const_cm","double",(void *)&_const_cm,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"\n\n\n\n// #line 108 \"cylinder.c\"\nstatic real dot3(const real *a, const real *b) {\n  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];\n}"," \n// #line 113 \"cylinder.c\"\n{\n    real delta;\n    delta = (2. * _const_cm * Delta + 0.);\n    real fx[3] = {val(fm.x,1,0,0) - val(fm.x,0,0,0), val(fm.x,1,0,0), -val(fm.x,0,0,0)};\n    real fy[3] = {val(fm.y,0, 1,0) - val(fm.y,0,0,0), val(fm.y,0, 1,0), -val(fm.y,0,0,0)};\n    real fz[3] = {val(fm.z,0, 0, 1) - val(fm.z,0,0,0), val(fm.z,0, 0, 1), -val(fm.z,0,0,0)};\n    real xy[3] = {val(u.x,0,0,0), val(u.x,0, 1,0), val(u.x,0, -1,0)};\n    real xz[3] = {val(u.x,0,0,0), val(u.x,0, 0, 1), val(u.x,0, 0, -1)};\n    real yx[3] = {val(u.y,0,0,0), val(u.y,1,0,0), val(u.y,-1,0,0)};\n    real yz[3] = {val(u.y,0,0,0), val(u.y,0, 0, 1), val(u.y,0, 0, -1)};\n    real zx[3] = {val(u.z,0,0,0), val(u.z,1,0,0), val(u.z,-1,0,0)};\n    real zy[3] = {val(u.z,0,0,0), val(u.z,0, 1,0), val(u.z,0, -1,0)};\n    val_out_(omega.x,0,0,0) = (dot3(fy, zy) - dot3(fz, yz)) / delta;\n    val_out_(omega.y,0,0,0) = (dot3(fz, xz) - dot3(fx, zx)) / delta;\n    val_out_(omega.z,0,0,0) = (dot3(fx, yx) - dot3(fy, xy)) / delta;\n  }"}) {

;      
          
    
#line 116
_stencil_val(fm.x,0,0,0); _stencil_val(fm.x,1,0,0); _stencil_val(fm.x,0,0,0);_stencil_val(fm.x,1,0,0);     
    _stencil_val(fm.y,0,0,0); _stencil_val(fm.y,0, 1,0); _stencil_val(fm.y,0,0,0);_stencil_val(fm.y,0, 1,0);     
    _stencil_val(fm.z,0,0,0); _stencil_val(fm.z,0, 0, 1); _stencil_val(fm.z,0,0,0);_stencil_val(fm.z,0, 0, 1);   
     _stencil_val(u.x,0, -1,0); _stencil_val(u.x,0, 1,0);_stencil_val(u.x,0,0,0);   
     _stencil_val(u.x,0, 0, -1); _stencil_val(u.x,0, 0, 1);_stencil_val(u.x,0,0,0);   
     _stencil_val(u.y,-1,0,0); _stencil_val(u.y,1,0,0);_stencil_val(u.y,0,0,0);   
     _stencil_val(u.y,0, 0, -1); _stencil_val(u.y,0, 0, 1);_stencil_val(u.y,0,0,0);   
     _stencil_val(u.z,-1,0,0); _stencil_val(u.z,1,0,0);_stencil_val(u.z,0,0,0);   
     _stencil_val(u.z,0, -1,0); _stencil_val(u.z,0, 1,0);_stencil_val(u.z,0,0,0);
    _stencil_val_a(omega.x,0,0,0);        
    _stencil_val_a(omega.y,0,0,0);        
    _stencil_val_a(omega.z,0,0,0);        
  }end_foreach_stencil();
  {
#line 113
foreach () {
    double delta;
    delta = (2. * _const_cm * Delta + 0.);
    double fx[3] = {val(fm.x,1,0,0) - val(fm.x,0,0,0), val(fm.x,1,0,0), -val(fm.x,0,0,0)};
    double fy[3] = {val(fm.y,0, 1,0) - val(fm.y,0,0,0), val(fm.y,0, 1,0), -val(fm.y,0,0,0)};
    double fz[3] = {val(fm.z,0, 0, 1) - val(fm.z,0,0,0), val(fm.z,0, 0, 1), -val(fm.z,0,0,0)};
    double xy[3] = {val(u.x,0,0,0), val(u.x,0, 1,0), val(u.x,0, -1,0)};
    double xz[3] = {val(u.x,0,0,0), val(u.x,0, 0, 1), val(u.x,0, 0, -1)};
    double yx[3] = {val(u.y,0,0,0), val(u.y,1,0,0), val(u.y,-1,0,0)};
    double yz[3] = {val(u.y,0,0,0), val(u.y,0, 0, 1), val(u.y,0, 0, -1)};
    double zx[3] = {val(u.z,0,0,0), val(u.z,1,0,0), val(u.z,-1,0,0)};
    double zy[3] = {val(u.z,0,0,0), val(u.z,0, 1,0), val(u.z,0, -1,0)};
    val(omega.x,0,0,0) = (dot3(fy, zy) - dot3(fz, yz)) / delta;
    val(omega.y,0,0,0) = (dot3(fz, xz) - dot3(fx, zx)) / delta;
    val(omega.z,0,0,0) = (dot3(fx, yx) - dot3(fy, xy)) / delta;
  }end_foreach();}}else if(!is_constant(cm) && is_constant(fm.x)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);
  
#line 113
foreach_stencil (1,{(NonLocal[]){{"omega","vector",(void *)&omega,NULL,0},{"u","vector",(void *)&u,NULL,0},{"_const_fm","_coord",(void *)&_const_fm,NULL,0},{"cm","scalar",(void *)&cm,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"\n\n\n\n// #line 108 \"cylinder.c\"\nstatic real dot3(const real *a, const real *b) {\n  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];\n}"," \n// #line 113 \"cylinder.c\"\n{\n    real delta;\n    delta = (2. * val(cm,0,0,0) * Delta + 0.);\n    real fx[3] = {_const_fm.x - _const_fm.x, _const_fm.x, -_const_fm.x};\n    real fy[3] = {_const_fm.y - _const_fm.y, _const_fm.y, -_const_fm.y};\n    real fz[3] = {_const_fm.z - _const_fm.z, _const_fm.z, -_const_fm.z};\n    real xy[3] = {val(u.x,0,0,0), val(u.x,0, 1,0), val(u.x,0, -1,0)};\n    real xz[3] = {val(u.x,0,0,0), val(u.x,0, 0, 1), val(u.x,0, 0, -1)};\n    real yx[3] = {val(u.y,0,0,0), val(u.y,1,0,0), val(u.y,-1,0,0)};\n    real yz[3] = {val(u.y,0,0,0), val(u.y,0, 0, 1), val(u.y,0, 0, -1)};\n    real zx[3] = {val(u.z,0,0,0), val(u.z,1,0,0), val(u.z,-1,0,0)};\n    real zy[3] = {val(u.z,0,0,0), val(u.z,0, 1,0), val(u.z,0, -1,0)};\n    val_out_(omega.x,0,0,0) = (dot3(fy, zy) - dot3(fz, yz)) / delta;\n    val_out_(omega.y,0,0,0) = (dot3(fz, xz) - dot3(fx, zx)) / delta;\n    val_out_(omega.z,0,0,0) = (dot3(fx, yx) - dot3(fy, xy)) / delta;\n  }"}) { 
     
_stencil_val(cm,0,0,0);
;;;;
;;;;
;;;;   
     _stencil_val(u.x,0, -1,0); _stencil_val(u.x,0, 1,0);_stencil_val(u.x,0,0,0);   
     _stencil_val(u.x,0, 0, -1); _stencil_val(u.x,0, 0, 1);_stencil_val(u.x,0,0,0);   
     _stencil_val(u.y,-1,0,0); _stencil_val(u.y,1,0,0);_stencil_val(u.y,0,0,0);   
     _stencil_val(u.y,0, 0, -1); _stencil_val(u.y,0, 0, 1);_stencil_val(u.y,0,0,0);   
     _stencil_val(u.z,-1,0,0); _stencil_val(u.z,1,0,0);_stencil_val(u.z,0,0,0);   
     _stencil_val(u.z,0, -1,0); _stencil_val(u.z,0, 1,0);_stencil_val(u.z,0,0,0);
    _stencil_val_a(omega.x,0,0,0);        
    _stencil_val_a(omega.y,0,0,0);        
    _stencil_val_a(omega.z,0,0,0);        
  }end_foreach_stencil();
  {
#line 113
foreach () {
    double delta;
    delta = (2. * val(cm,0,0,0) * Delta + 0.);
    double fx[3] = {_const_fm.x - _const_fm.x, _const_fm.x, -_const_fm.x};
    double fy[3] = {_const_fm.y - _const_fm.y, _const_fm.y, -_const_fm.y};
    double fz[3] = {_const_fm.z - _const_fm.z, _const_fm.z, -_const_fm.z};
    double xy[3] = {val(u.x,0,0,0), val(u.x,0, 1,0), val(u.x,0, -1,0)};
    double xz[3] = {val(u.x,0,0,0), val(u.x,0, 0, 1), val(u.x,0, 0, -1)};
    double yx[3] = {val(u.y,0,0,0), val(u.y,1,0,0), val(u.y,-1,0,0)};
    double yz[3] = {val(u.y,0,0,0), val(u.y,0, 0, 1), val(u.y,0, 0, -1)};
    double zx[3] = {val(u.z,0,0,0), val(u.z,1,0,0), val(u.z,-1,0,0)};
    double zy[3] = {val(u.z,0,0,0), val(u.z,0, 1,0), val(u.z,0, -1,0)};
    val(omega.x,0,0,0) = (dot3(fy, zy) - dot3(fz, yz)) / delta;
    val(omega.y,0,0,0) = (dot3(fz, xz) - dot3(fx, zx)) / delta;
    val(omega.z,0,0,0) = (dot3(fx, yx) - dot3(fy, xy)) / delta;
  }end_foreach();}}else {double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);
  
#line 113
foreach_stencil (1,{(NonLocal[]){{"omega","vector",(void *)&omega,NULL,0},{"u","vector",(void *)&u,NULL,0},{"_const_fm","_coord",(void *)&_const_fm,NULL,0},{"_const_cm","double",(void *)&_const_cm,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"\n\n\n\n// #line 108 \"cylinder.c\"\nstatic real dot3(const real *a, const real *b) {\n  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];\n}"," \n// #line 113 \"cylinder.c\"\n{\n    real delta;\n    delta = (2. * _const_cm * Delta + 0.);\n    real fx[3] = {_const_fm.x - _const_fm.x, _const_fm.x, -_const_fm.x};\n    real fy[3] = {_const_fm.y - _const_fm.y, _const_fm.y, -_const_fm.y};\n    real fz[3] = {_const_fm.z - _const_fm.z, _const_fm.z, -_const_fm.z};\n    real xy[3] = {val(u.x,0,0,0), val(u.x,0, 1,0), val(u.x,0, -1,0)};\n    real xz[3] = {val(u.x,0,0,0), val(u.x,0, 0, 1), val(u.x,0, 0, -1)};\n    real yx[3] = {val(u.y,0,0,0), val(u.y,1,0,0), val(u.y,-1,0,0)};\n    real yz[3] = {val(u.y,0,0,0), val(u.y,0, 0, 1), val(u.y,0, 0, -1)};\n    real zx[3] = {val(u.z,0,0,0), val(u.z,1,0,0), val(u.z,-1,0,0)};\n    real zy[3] = {val(u.z,0,0,0), val(u.z,0, 1,0), val(u.z,0, -1,0)};\n    val_out_(omega.x,0,0,0) = (dot3(fy, zy) - dot3(fz, yz)) / delta;\n    val_out_(omega.y,0,0,0) = (dot3(fz, xz) - dot3(fx, zx)) / delta;\n    val_out_(omega.z,0,0,0) = (dot3(fx, yx) - dot3(fy, xy)) / delta;\n  }"}) {

;
;;;;
;;;;
;;;;   
     _stencil_val(u.x,0, -1,0); _stencil_val(u.x,0, 1,0);_stencil_val(u.x,0,0,0);   
     _stencil_val(u.x,0, 0, -1); _stencil_val(u.x,0, 0, 1);_stencil_val(u.x,0,0,0);   
     _stencil_val(u.y,-1,0,0); _stencil_val(u.y,1,0,0);_stencil_val(u.y,0,0,0);   
     _stencil_val(u.y,0, 0, -1); _stencil_val(u.y,0, 0, 1);_stencil_val(u.y,0,0,0);   
     _stencil_val(u.z,-1,0,0); _stencil_val(u.z,1,0,0);_stencil_val(u.z,0,0,0);   
     _stencil_val(u.z,0, -1,0); _stencil_val(u.z,0, 1,0);_stencil_val(u.z,0,0,0);
    _stencil_val_a(omega.x,0,0,0);        
    _stencil_val_a(omega.y,0,0,0);        
    _stencil_val_a(omega.z,0,0,0);        
  }end_foreach_stencil();
  {
#line 113
foreach () {
    double delta;
    delta = (2. * _const_cm * Delta + 0.);
    double fx[3] = {_const_fm.x - _const_fm.x, _const_fm.x, -_const_fm.x};
    double fy[3] = {_const_fm.y - _const_fm.y, _const_fm.y, -_const_fm.y};
    double fz[3] = {_const_fm.z - _const_fm.z, _const_fm.z, -_const_fm.z};
    double xy[3] = {val(u.x,0,0,0), val(u.x,0, 1,0), val(u.x,0, -1,0)};
    double xz[3] = {val(u.x,0,0,0), val(u.x,0, 0, 1), val(u.x,0, 0, -1)};
    double yx[3] = {val(u.y,0,0,0), val(u.y,1,0,0), val(u.y,-1,0,0)};
    double yz[3] = {val(u.y,0,0,0), val(u.y,0, 0, 1), val(u.y,0, 0, -1)};
    double zx[3] = {val(u.z,0,0,0), val(u.z,1,0,0), val(u.z,-1,0,0)};
    double zy[3] = {val(u.z,0,0,0), val(u.z,0, 1,0), val(u.z,0, -1,0)};
    val(omega.x,0,0,0) = (dot3(fy, zy) - dot3(fz, yz)) / delta;
    val(omega.y,0,0,0) = (dot3(fz, xz) - dot3(fx, zx)) / delta;
    val(omega.z,0,0,0) = (dot3(fx, yx) - dot3(fy, xy)) / delta;
  }end_foreach();}}
}

static int slice(double x, double y, double z, double Delta) {
  double epsilon = Delta / 10;
  return z <= -epsilon && z + Delta + epsilon >= 0;
}
static double shape_cylinder(double x, double y, double z) {
  return sq(x) + sq(y) - sq(1.0 / 2);
}
static double shape_sphere(double x, double y, double z) {
  return sq(x) + sq(y) + sq(z) - sq(1.0 / 2);
}
static double (*Shape[])(double, double, double) = {shape_cylinder,
                                                    shape_sphere};
static const char *shape_names[] = {"cylinder", "sphere"};
static double (*shape)(double, double, double);

static const char *force_path, *output_prefix;
static char *dump_path;
static const int outlevel = 5;
static double reynolds, tend;
static int maxlevel, minlevel, period, Verbose, FullOutput, AdaptFlag;
static vector  muv={{15},{16},{17}};
static scalar  l2={18};
static vector  omega={{19},{20},{21}};
static scalar  phi={22};

static double _boundary6(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _dirichlet(1, point, neighbor, _s, data));}}}static double _boundary6_homogeneous(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _dirichlet_homogeneous(1, point, neighbor, _s, data));}}}
static double _boundary7(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann(0, point, neighbor, _s, data));}}}static double _boundary7_homogeneous(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann_homogeneous(0, point, neighbor, _s, data));}}}
static double _boundary8(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann(0, point, neighbor, _s, data));}}}static double _boundary8_homogeneous(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann_homogeneous(0, point, neighbor, _s, data));}}}

static double _boundary9(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann(0, point, neighbor, _s, data));}}}static double _boundary9_homogeneous(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann_homogeneous(0, point, neighbor, _s, data));}}}
static double _boundary10(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _dirichlet(0, point, neighbor, _s, data));}}}static double _boundary10_homogeneous(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _dirichlet_homogeneous(0, point, neighbor, _s, data));}}}
static double _boundary11(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _dirichlet(0, point, neighbor, _s, data));}}}static double _boundary11_homogeneous(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _dirichlet_homogeneous(0, point, neighbor, _s, data));}}}

int main(int argc, char **argv) {
#line 443
_init_solver();
  
#line 165
char *end;
  const char *periodic_boundaries;
  int ReynoldsFlag, MaxLevelFlag, MinLevelFlag, PeriodFlag, TendFlag,
      DomainFlag, DTFlag, i;
  double domain, dt_min;
  AdaptFlag = 0;
  DomainFlag = 0;
  FullOutput = 0;
  MaxLevelFlag = 0;
  MinLevelFlag = 0;
  PeriodFlag = 0;
  ReynoldsFlag = 0;
  TendFlag = 0;
  DTFlag = 0;
  Verbose = 0;
  output_prefix = NULL;
  force_path = NULL;
  dump_path = NULL;
  shape = NULL;
  periodic_boundaries = NULL;
  while (*++argv != NULL && argv[0][0] == '-')
    switch (argv[0][1]) {
    case 'h':
      fprintf(
          ferr,
          "Usage: cylinder [-h] [-v] [-a] [-F] [-r <Reynolds number> "
          "-l <resolution level> -m <maximum resolution level> "
          "-o <prefix> -p <dump period> -e <end time> "
          "-f <force file> -S cylinder|sphere "
          "-z <domain size> [-b <boundaries>] [-d <dump file>]\n\n"
          "Options:\n"
          "  -h          Display this help message\n"
          "  -v          Verbose\n"
          "  -F          Output the full field\n"
          "  -a          Use adoptation\n"
          "  -b <string> Periodic boundary (ft|f|t: front (f), top (t) or both,"
          "default is symmetric boundary)\n"
          "  -r <num>    Reynolds number\n"
          "  -t <num>    dt\n"
          "  -l <num>    Minimum resolution level (positive integer)\n"
          "  -m <num>    Maximum resolution level (positive integer)\n"
          "  -o <string> Prefix for the output files\n"
          "  -p <num>    Dump period (positive integer)\n"
          "  -e <num>    End time of the simulation (decimal number)\n"
          "  -f <file>   Output force file\n"
          "  -S <string> Specify shape (cylinder|sphere)\n"
          "  -d <file>   Restart simulation from the dump file\n"
          "  -z <num>    Domain size\n\n"
   "  npe: %d\n"
          "Example usage:\n"
          "  ./cylinder -v -r 100 -l 7 -m 10 -p 100 -e 2 -z 2.5 -S sphere\n"
          "  ./cylinder -v -r 100 -l 7 -m 10 -p 100 -e 2 -f force.dat -z 2.5 "
          "-S cylinder -o h -b t\n", npe());
      exit(1);
    case 'r':
      argv++;
      if (*argv == NULL) {
        fprintf(ferr, "cylinder: error: -r needs an argument\n");
        exit(1);
      }
      reynolds = strtod(*argv, &end);
      if (*end != '\0') {
        fprintf(ferr, "cylinder: error: '%s' is not a number\n", *argv);
        exit(1);
      }
      ReynoldsFlag = 1;
      break;
    case 't':
      argv++;
      if (*argv == NULL) {
        fprintf(ferr, "cylinder: error: -t needs an argument\n");
        exit(1);
      }
      dt_min = strtod(*argv, &end);
      if (*end != '\0') {
        fprintf(ferr, "cylinder: error: '%s' is not a number\n", *argv);
        exit(1);
      }
      DTFlag = 1;
      break;
    case 'm':
      argv++;
      if (*argv == NULL) {
        fprintf(ferr, "cylinder: error: -m needs an argument\n");
        exit(1);
      }
      maxlevel = strtol(*argv, &end, 10);
      if (*end != '\0' || maxlevel <= 0) {
        fprintf(ferr, "cylinder: error: '%s' is not a positive integer\n",
                *argv);
        exit(1);
      }
      MaxLevelFlag = 1;
      break;
    case 'l':
      argv++;
      if (*argv == NULL) {
        fprintf(ferr, "cylinder: error: -l needs an argument\n");
        exit(1);
      }
      minlevel = strtol(*argv, &end, 10);
      if (*end != '\0' || minlevel <= 0) {
        fprintf(ferr, "cylinder: error: '%s' is not a positive integer\n",
                *argv);
        exit(1);
      }
      MinLevelFlag = 1;
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
    case 'v':
      Verbose = 1;
      break;
    case 'a':
      AdaptFlag = 1;
      break;
    case 'd':
      argv++;
      if (*argv == NULL) {
        fprintf(ferr, "cylinder: error: -d needs an argument\n");
        exit(1);
      }
      dump_path = *argv;
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
    case 'f':
      argv++;
      if (*argv == NULL) {
        fprintf(ferr, "cylinder: error: -f needs an argument\n");
        exit(1);
      }
      force_path = *argv;
      break;
    case 'F':
      FullOutput = 1;
      break;
    case 'S':
      argv++;
      if (*argv == NULL) {
        fprintf(ferr, "cylinder: error: -S needs an argument\n");
        exit(1);
      }
      for (i = 0; ; i++) {
        if (i == sizeof shape_names / sizeof *shape_names) {
          fprintf(ferr, "cylinder: error: unknown shape '%s'\n", *argv);
          exit(1);
        }
        if (strcmp(shape_names[i], *argv) == 0) {
          shape = Shape[i];
          break;
        }
      }
      break;
    case 'b':
      argv++;
      if (*argv == NULL) {
        fprintf(ferr, "cylinder: error: -b needs an argument\n");
        exit(1);
      }
      periodic_boundaries = *argv;
      break;
    case 'o':
      argv++;
      if (*argv == NULL) {
        fprintf(ferr, "cylinder: error: -o needs an argument\n");
        exit(1);
      }
      output_prefix = *argv;
      break;
    case 'z':
      argv++;
      if (*argv == NULL) {
        fprintf(ferr, "cylinder: error: -z needs an argument\n");
        exit(1);
      }
      domain = strtod(*argv, &end);
      if (*end != '\0') {
        fprintf(ferr, "cylinder: error: '%s' is not a number\n", *argv);
        exit(1);
      }
      if (domain < 1) {
        fprintf(ferr,
                "cylinder: error: '%s': domain size (-z) is less then one\n",
                *argv);
        exit(1);
      }
      DomainFlag = 1;
      break;
    default:
      fprintf(ferr, "cylinder: error: unknown option '%s'\n", *argv);
      exit(1);
    }
  if (!ReynoldsFlag) {
    fprintf(ferr, "cylinder: error: -r is not set\n");
    exit(1);
  }
  if (!MinLevelFlag) {
    fprintf(ferr, "cylinder: error: -l is not set\n");
    exit(1);
  }
  if (!MaxLevelFlag) {
    fprintf(ferr, "cylinder: error: -m is not set\n");
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
  if (dump_path == NULL && !DomainFlag) {
    fprintf(ferr,
            "cylinder: error: eather -d (dump) -z (size) must be set\n");
    exit(1);
  }
  if (dump_path == NULL && shape == NULL) {
    fprintf(ferr,
            "cylinder: error: eather -d (dump) or -S (shape) must be set\n");
    exit(1);
  }
  if (Verbose && pid() == 0)
    fprintf(ferr, "cylinder: starting on %d ranks\n", npe());
  if (dump_path == NULL) {
    size(domain);
    origin(-L0 / 2.5, -L0 / 2.0, -L0 / 2.0);
  }
  mu = muv;
  if (periodic_boundaries != NULL)
    for (i = 0; periodic_boundaries[i] != '\0'; i++)
      switch (periodic_boundaries[i]) {
      case 'f':
        tree_periodic(front);
        if (Verbose && pid() == 0)
          fprintf(ferr, "cylinder: front boundary is periodic\n");
        break;
      case 't':
        tree_periodic(top);
        if (Verbose && pid() == 0)
          fprintf(ferr, "cylinder: top boundary is periodic\n");
        break;
      default:
        fprintf(ferr, "cylinder: unknown boundary in '%s'\n",
                periodic_boundaries);
        exit(1);
        break;
      }
  if (DTFlag)
    DT = dt_min;
  init_grid(1 << outlevel);
  run();
  if (Verbose && pid() == 0)
    fprintf(ferr, "cylinder: done\n");
free_solver();

#line 443
}

static int init_0_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(t = 0)!=0;*ip=i;*tp=t;return ret;}


#line 445
      static int init_0(const int i,const double t,Event *_ev){tracing("init_0","cylinder.c",445); {
  FILE *dump_file;
  if (dump_path == NULL) {
    do { int refined; do { boundary_internal ((scalar *)all, "cylinder.c", 448); refined = 0; ((Tree *)grid)->refined.n = 0; {foreach_leaf() if (x < X0 + 0.9 * L0 && level < minlevel) { refine_cell (point, all, 0, &((Tree *)grid)->refined); refined++; continue; }end_foreach_leaf();} mpi_all_reduce (refined, MPI_INT, MPI_SUM); if (refined) { mpi_boundary_refine (all); mpi_boundary_update (all); } } while (refined); } while(0);
    for (;;) {
      do { scalar  phi=new_vertex_scalar("phi"); foreach_vertex_stencil(1,{(NonLocal[]){{"phi","scalar",(void *)&phi,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},""," \n// #line 450 \"cylinder.c\"\nval_out_(phi,0,0,0) = shape(x, y, z);"}) {_stencil_val_a(phi,0,0,0);    }end_foreach_vertex_stencil(); {foreach_vertex() val(phi,0,0,0) = shape(x, y, z);end_foreach_vertex();} fractions (phi, cs, fs
#line 122 "/home/lisergey/basilisk/src/fractions.h"
, 0.
#line 450 "cylinder.c"
);delete((scalar*)((scalar[]){phi,{-1}})); } while(0);
      astats s = adapt_wavelet(
#line 451 "/home/lisergey/basilisk/src/grid/tree-common.h"
(
#line 173
scalar *
#line 451
)
#line 451 "cylinder.c"
((scalar[]){cs,{-1}}), (double[]){0}, maxlevel
, minlevel
#line 176 "/home/lisergey/basilisk/src/grid/tree-common.h"
, 
all
#line 452 "cylinder.c"
);
      if (Verbose && pid() == 0)
        fprintf(ferr, "cylinder: refined %d cells\n", s.nf);
      if (s.nf == 0)
        break;
    }
  } else {
    if (Verbose && pid() == 0)
      fprintf(ferr, "cylinder: reading dump from '%s'\n", dump_path);
    if ((dump_file = fopen(dump_path, "r")) == NULL) {
      fprintf(ferr, "cylinder: error: failed to open '%s'\n", dump_path);
      exit(1);
    }
    restore( 
#line 1116 "/home/lisergey/basilisk/src/output.h"
"dump", 
NULL
#line 465 "cylinder.c"
, dump_file);





    if (Verbose && pid() == 0)
      fprintf(ferr, "cylinder: starting from '%s': time: %g, step: %d\n",
              dump_path, t, i);
    if (i == 0)
      fractions(phi, cs, fs
#line 122 "/home/lisergey/basilisk/src/fractions.h"
, 0.
#line 475 "cylinder.c"
);
    fractions_cleanup((struct Cleanup){cs, fs});
    if (Verbose)
      fields_stats();
  }
  if (i == 0) {
    if (Verbose && pid() == 0)
      fprintf(ferr, "cylinder: initialize velocity\n");
    foreach_stencil (1,{(NonLocal[]){{"cs","scalar",(void *)&cs,NULL,0},{"u","vector",(void *)&u,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},""," \n// #line 483 \"cylinder.c\"\n{\n      val_out_(u.x,0,0,0) = val(cs,0,0,0);\n      val_out_(u.y,0,0,0) = 0;\n      val_out_(u.z,0,0,0) = 0;\n    }"}) { 
_stencil_val(cs,0,0,0);
      
#line 484
_stencil_val_a(u.x,0,0,0); 
      _stencil_val_a(u.y,0,0,0);  
      _stencil_val_a(u.z,0,0,0);  
    }end_foreach_stencil();
    {
#line 483
foreach () {
      val(u.x,0,0,0) = val(cs,0,0,0);
      val(u.y,0,0,0) = 0;
      val(u.z,0,0,0) = 0;
    }end_foreach();}

    event("dump");
  }
}{end_tracing("init_0","cylinder.c",491);return 0;}end_tracing("init_0","cylinder.c",491);}

static int properties_0_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i++)!=0;*ip=i;*tp=t;return ret;}


#line 493
      static int properties_0(const int i,const double t,Event *_ev){tracing("properties_0","cylinder.c",493); { if(!is_constant(fm.x)){ foreach_face_stencil(1,{(NonLocal[]){{"reynolds","double",(void *)&reynolds,NULL,0},{"fm","vector",(void *)&fm,NULL,0},{"muv","vector",(void *)&muv,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n// #line 493 \"cylinder.c\"\n{is_face_x(){ val_out_(muv.x,0,0,0) = val(fm.x,0,0,0) / reynolds;}end_is_face_x()is_face_y(){ val_out_(muv.y,0,0,0) = val(fm.y,0,0,0) / reynolds;}end_is_face_y()is_face_z(){ val_out_(muv.z,0,0,0) = val(fm.z,0,0,0) / reynolds;}end_is_face_z()}"}){_stencil_is_face_x(){ { _stencil_val(fm.x,0,0,0);_stencil_val_a(muv.x,0,0,0);   }}end__stencil_is_face_x()_stencil_is_face_y(){ { _stencil_val(fm.y,0,0,0);_stencil_val_a(muv.y,0,0,0);   }}end__stencil_is_face_y()_stencil_is_face_z(){ { _stencil_val(fm.z,0,0,0);_stencil_val_a(muv.z,0,0,0);   }}end__stencil_is_face_z()}end_foreach_face_stencil();{foreach_face_generic(){is_face_x(){ val(muv.x,0,0,0) = val(fm.x,0,0,0) / reynolds;}end_is_face_x()is_face_y(){ val(muv.y,0,0,0) = val(fm.y,0,0,0) / reynolds;}end_is_face_y()is_face_z(){ val(muv.z,0,0,0) = val(fm.z,0,0,0) / reynolds;}end_is_face_z()}end_foreach_face_generic();}}else {_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm); foreach_face_stencil(1,{(NonLocal[]){{"reynolds","double",(void *)&reynolds,NULL,0},{"_const_fm","_coord",(void *)&_const_fm,NULL,0},{"muv","vector",(void *)&muv,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},"","\n// #line 493 \"cylinder.c\"\n{is_face_x(){ val_out_(muv.x,0,0,0) = _const_fm.x / reynolds;}end_is_face_x()is_face_y(){ val_out_(muv.y,0,0,0) = _const_fm.y / reynolds;}end_is_face_y()is_face_z(){ val_out_(muv.z,0,0,0) = _const_fm.z / reynolds;}end_is_face_z()}"}){_stencil_is_face_x(){ {;_stencil_val_a(muv.x,0,0,0);   }}end__stencil_is_face_x()_stencil_is_face_y(){ {;_stencil_val_a(muv.y,0,0,0);   }}end__stencil_is_face_y()_stencil_is_face_z(){ {;_stencil_val_a(muv.z,0,0,0);   }}end__stencil_is_face_z()}end_foreach_face_stencil(); {foreach_face_generic(){is_face_x(){ val(muv.x,0,0,0) = _const_fm.x / reynolds;}end_is_face_x()is_face_y(){ val(muv.y,0,0,0) = _const_fm.y / reynolds;}end_is_face_y()is_face_z(){ val(muv.z,0,0,0) = _const_fm.z / reynolds;}end_is_face_z()}end_foreach_face_generic();}} }{end_tracing("properties_0","cylinder.c",493);return 0;}end_tracing("properties_0","cylinder.c",493);}

static int dump_0_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=( t <= tend)!=0;*ip=i;*tp=t;return ret;}static int dump_0_expr1(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i++)!=0;*ip=i;*tp=t;return ret;}


#line 495
      static int dump_0(const int i,const double t,Event *_ev){tracing("dump_0","cylinder.c",495); {
  char path[FILENAME_MAX];
  coord Fp, Fmu;
  static FILE *fp;

  if (i % period == 0) {
    if (Verbose) {
      fields_stats();
      if (pid() == 0)
        fprintf(ferr, "cylinder: %d: %09d %.16e %ld\n", npe(), i, t, grid->n);
    }
    if (output_prefix != NULL) {
      vorticity_vector(u, omega);
      lambda2(u, l2);
      if (FullOutput) {
        sprintf(path, "%s.%09d", output_prefix, i);
        output_xdmf(t,((scalar[]) {p, l2,{-1}}),((vector[]) {u, omega,{{-1},{-1},{-1}}}), NULL, path);
      }
      snprintf(path, sizeof path, "%s.slice.%09d", output_prefix, i);
      output_xdmf(t,((scalar[]) {p, l2, cs, phi,{-1}}),((vector[]) {u, omega,{{-1},{-1},{-1}}}), slice, path);
      if (i % (10 * period) == 0) {
        snprintf(path, sizeof path, "%s.%09d.dump", output_prefix, i);
        dump(path
#line 1041 "/home/lisergey/basilisk/src/output.h"
, 
all, 
NULL, 
false
#line 517 "cylinder.c"
);
      }
    }
    if (force_path) {
      embed_force3(p, u, mu, &Fp, &Fmu);
      if (pid() == 0) {
        if (fp == NULL) {
          if ((fp = fopen(force_path, "w")) == NULL) {
            fprintf(ferr, "cylinder: error: fail to open '%s'\n", force_path);
            exit(1);
          }
        } else {
          if ((fp = fopen(force_path, "a")) == NULL) {
            fprintf(ferr, "cylinder: error: fail to open '%s'\n", force_path);
            exit(1);
          }
        }
        fprintf(fp,
                "%d %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e "
                "%.16e %.16e\n",
                i, t, Fp.x + Fmu.x, Fp.y + Fmu.y, Fp.z + Fmu.z, Fp.x, Fp.y,
                Fp.z, Fmu.x, Fmu.y, Fmu.z, dt);
        fflush(fp);
      }
    }
  }
  if (AdaptFlag) {
    astats s = adapt_wavelet((scalar *)((scalar[]){cs, u.x, u.y, u.z,{-1}}), (double[]){0, 0.1, 0.1, 0.1}
, maxlevel, minlevel
#line 176 "/home/lisergey/basilisk/src/grid/tree-common.h"
, 
all
#line 545 "cylinder.c"
);

    do { static const int too_fine = 1 << user; {foreach_cell() { if (is_leaf(cell)) continue; if (is_local(cell) && (!(x < X0 + 0.9 * L0) && level > outlevel)) cell.flags |= too_fine; }end_foreach_cell();} for (int _l = depth(); _l >= 0; _l--) { {foreach_cell() { if (is_leaf(cell)) continue; if (level == _l) { if (is_local(cell) && (cell.flags & too_fine)) { coarsen_cell (point, all); cell.flags &= ~too_fine; } continue; } }end_foreach_cell();} mpi_boundary_coarsen (_l, too_fine); } mpi_boundary_update (all); } while (0);
    fractions_cleanup((struct Cleanup){cs, fs});
    if (Verbose && i % period == 0 && pid() == 0)
      fprintf(ferr, "cylinder: refined %d cells, coarsened %d cells\n", s.nf,
              s.nc);
  }
}{end_tracing("dump_0","cylinder.c",553);return 0;}end_tracing("dump_0","cylinder.c",553);}
#line 2 "ast/init_solver.h"

static void _init_solver (void)
{
  void init_solver();
datasize=23*sizeof(double);
  
#line 6
init_solver();
  {
#line 24
octree_methods();

    

    
#line 12
{
      
      
    
      {  
#line 42 "/home/lisergey/basilisk/src/run.h"
event_register((Event){0,1,defaults,{defaults_expr0},((int *)0),((double *)0),"/home/lisergey/basilisk/src/run.h",42,"defaults"});  
#line 126 "/home/lisergey/basilisk/src/navier-stokes/centered.h"
event_register((Event){0,1,defaults_0,{defaults_0_expr0},((int *)0),((double *)0),"/home/lisergey/basilisk/src/navier-stokes/centered.h",126,"defaults"});  
#line 187
event_register((Event){0,1,default_display,{default_display_expr0},((int *)0),((double *)0),"/home/lisergey/basilisk/src/navier-stokes/centered.h",187,"default_display"});  








event_register((Event){0,1,init,{init_expr0},((int *)0),((double *)0),"/home/lisergey/basilisk/src/navier-stokes/centered.h",196,"init"});  
#line 96 "cylinder.c"
event_register((Event){0,1,velocity,{velocity_expr0},((int *)0),((double *)0),"cylinder.c",96,"velocity"});  
#line 445
event_register((Event){0,1,init_0,{init_0_expr0},((int *)0),((double *)0),"cylinder.c",445,"init"});  
#line 495
event_register((Event){0,2,dump_0,{dump_0_expr0,dump_0_expr1},((int *)0),((double *)0),"cylinder.c",495,"dump"});
	
	
	
      
#line 22 "ast/init_solver.h"
}
#line 1257 "/home/lisergey/basilisk/src/common.h"
init_const_vector((vector){{_NVARMAX+0},{_NVARMAX+1},{_NVARMAX+2}},"zerof",(double[]){0.,0.,0.});
init_const_vector((vector){{_NVARMAX+3},{_NVARMAX+4},{_NVARMAX+5}},"unityf",(double[]){1.,1.,1.});
init_const_scalar((scalar){_NVARMAX+6},"unity", 1.);
init_const_scalar((scalar){_NVARMAX+7},"zeroc", 0.);



init_const_vector((vector){{_NVARMAX+8},{_NVARMAX+9},{_NVARMAX+10}},"unityf0",(double[]){1.,1.,1.});
init_const_scalar((scalar){_NVARMAX+11},"unity0", 1.);  init_scalar((scalar){0},"p");  init_vector((vector){{1},{2},{3}},"u");  init_vector((vector){{4},{5},{6}},"g");  init_scalar((scalar){7},"pf");  init_face_vector((vector){{8},{9},{10}},"uf");  init_scalar((scalar){11},"cs");  init_face_vector((vector){{12},{13},{14}},"fs");  init_face_vector((vector){{15},{16},{17}},"muv");  init_scalar((scalar){18},"l2");  init_vector((vector){{19},{20},{21}},"omega");  init_scalar((scalar){22},"phi");
    
#line 23 "ast/init_solver.h"
}_attribute[p.i].dirty=1,_attribute[p.i].boundary[right]=_boundary0,_attribute[p.i].boundary_homogeneous[right]=_boundary0_homogeneous;_attribute[p.i].dirty=1,_attribute[p.i].boundary[left]=_boundary1,_attribute[p.i].boundary_homogeneous[left]=_boundary1_homogeneous;_attribute[p.i].dirty=1,_attribute[p.i].boundary[top]=_boundary2,_attribute[p.i].boundary_homogeneous[top]=_boundary2_homogeneous;_attribute[p.i].dirty=1,_attribute[p.i].boundary[bottom]=_boundary3,_attribute[p.i].boundary_homogeneous[bottom]=_boundary3_homogeneous;_attribute[p.i].dirty=1,_attribute[p.i].boundary[front]=_boundary4,_attribute[p.i].boundary_homogeneous[front]=_boundary4_homogeneous;_attribute[p.i].dirty=1,_attribute[p.i].boundary[back]=_boundary5,_attribute[p.i].boundary_homogeneous[back]=_boundary5_homogeneous;_attribute[u.x.i].dirty=1,_attribute[u.x.i].boundary[left]=_boundary6,_attribute[u.x.i].boundary_homogeneous[left]=_boundary6_homogeneous;_attribute[p.i].dirty=1,_attribute[p.i].boundary[left]=_boundary7,_attribute[p.i].boundary_homogeneous[left]=_boundary7_homogeneous;_attribute[pf.i].dirty=1,_attribute[pf.i].boundary[left]=_boundary8,_attribute[pf.i].boundary_homogeneous[left]=_boundary8_homogeneous;_attribute[u.x.i].dirty=1,_attribute[u.x.i].boundary[right]=_boundary9,_attribute[u.x.i].boundary_homogeneous[right]=_boundary9_homogeneous;_attribute[p.i].dirty=1,_attribute[p.i].boundary[right]=_boundary10,_attribute[p.i].boundary_homogeneous[right]=_boundary10_homogeneous;_attribute[pf.i].dirty=1,_attribute[pf.i].boundary[right]=_boundary11,_attribute[pf.i].boundary_homogeneous[right]=_boundary11_homogeneous;  
#line 50 "/home/lisergey/basilisk/src/run.h"
event_register((Event){0,1,cleanup,{cleanup_expr0},((int *)0),((double *)0),"/home/lisergey/basilisk/src/run.h",50,"cleanup"});  
#line 222 "/home/lisergey/basilisk/src/navier-stokes/centered.h"
event_register((Event){0,1,set_dtmax,{set_dtmax_expr0},((int *)0),((double *)0),"/home/lisergey/basilisk/src/navier-stokes/centered.h",222,"set_dtmax"});  

event_register((Event){0,1,stability,{stability_expr0},((int *)0),((double *)0),"/home/lisergey/basilisk/src/navier-stokes/centered.h",224,"stability"});  









event_register((Event){0,1,vof,{vof_expr0},((int *)0),((double *)0),"/home/lisergey/basilisk/src/navier-stokes/centered.h",234,"vof"});  
event_register((Event){0,1,tracer_advection,{tracer_advection_expr0},((int *)0),((double *)0),"/home/lisergey/basilisk/src/navier-stokes/centered.h",235,"tracer_advection"});  
event_register((Event){0,1,tracer_diffusion,{tracer_diffusion_expr0},((int *)0),((double *)0),"/home/lisergey/basilisk/src/navier-stokes/centered.h",236,"tracer_diffusion"});  






event_register((Event){0,1,properties,{properties_expr0},((int *)0),((double *)0),"/home/lisergey/basilisk/src/navier-stokes/centered.h",243,"properties"});  
#line 316
event_register((Event){0,1,advection_term,{advection_term_expr0},((int *)0),((double *)0),"/home/lisergey/basilisk/src/navier-stokes/centered.h",316,"advection_term"});  
#line 345
event_register((Event){0,1,viscous_term,{viscous_term_expr0},((int *)0),((double *)0),"/home/lisergey/basilisk/src/navier-stokes/centered.h",345,"viscous_term"});  
#line 381
event_register((Event){0,1,acceleration,{acceleration_expr0},((int *)0),((double *)0),"/home/lisergey/basilisk/src/navier-stokes/centered.h",381,"acceleration"});  
#line 421
event_register((Event){0,1,projection,{projection_expr0},((int *)0),((double *)0),"/home/lisergey/basilisk/src/navier-stokes/centered.h",421,"projection"});  
#line 436
event_register((Event){0,1,end_timestep,{end_timestep_expr0},((int *)0),((double *)0),"/home/lisergey/basilisk/src/navier-stokes/centered.h",436,"end_timestep"});  









event_register((Event){0,1,adapt,{adapt_expr0},((int *)0),((double *)0),"/home/lisergey/basilisk/src/navier-stokes/centered.h",446,"adapt"});  
#line 493 "cylinder.c"
event_register((Event){0,1,properties_0,{properties_0_expr0},((int *)0),((double *)0),"cylinder.c",493,"properties"});
  
#line 24 "ast/init_solver.h"
}
  set_fpe();
}
