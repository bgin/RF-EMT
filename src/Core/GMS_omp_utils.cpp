
#include <cstdio>
#include <cstdlib>
#include <sched.h>
#include <cerrno>
#include <sys/types.h>
#include <sys/resource.h> //getpid 
#include <signal.h>
#include <unistd.h>
#include <immintrin.h>
#include <omp.h> 
#include "GMS_omp_utils.h"

std::int32_t 
gms::common
::set_thread_affinity_and_priority(const std::int32_t cpu, const std::int32_t priority)
{
    cpu_set_t cpu_set;
    sched_param sp;
    int32_t status{-1};
    CPU_ZERO(&cpu_set);
    CPU_SET(cpu,&cpu_set);
    if(sched_setaffinity(0,sizeof(cpu_set), &cpu_set) < 0) 
    {
         status = 1;
         return status;
    }
    [[maybe_unused]] std::int32_t printf_ret = std::printf("Affinity set to cpu: %d\n",cpu);
    __builtin_memset(&sp,0,sizeof(sp));
    sp.sched_priority = priority; //99
    if((sched_setscheduler(0,SCHED_FIFO,&sp)) == -1) 
    {
        status = 2;
        return status;
    }
    status = 0;
    return status;
}

void 
gms::common
::print_omp_thread_affinity()
{
#define BUFFER_STORE_SIZE 80
#define FORMAT_STORE_SIZE 80
    char default_format[FORMAT_STORE_SIZE];
    char format_specifier[] = "host=%20H tid=%0.4n binds_to=%A";
    char buffer[BUFFER_STORE_SIZE];
    std::size_t nchars{};
    std::size_t diff;
    nchars = omp_get_affinity_format(default_format,(std::size_t)FORMAT_STORE_SIZE);
    diff   = nchars-(std::size_t)FORMAT_STORE_SIZE;
    if(diff>0ull)
        nchars += diff;
    omp_set_affinity_format(format_specifier);
    nchars = omp_capture_affinity(&buffer[0],(std::size_t)BUFFER_STORE_SIZE,NULL);
    [[maybe_unused]] std::int32_t printf_ret = printf("tid=%d affinity:%s\n",omp_get_thread_num(),buffer);
    
}

/*
   Usage: 
   int cl_thread_limit;

  if (argc != 2) {
    usage();
    return 1;
  }
  cl_thread_limit = atoi(argv[1]);

  omp_set_dynamic(0);
  if (omp_get_thread_limit() != cl_thread_limit) {
    fprintf(stderr, "omp_get_thread_limit failed with %d, should be%d\n",
            omp_get_thread_limit(), cl_thread_limit);
    return 1;
  }
  else if (omp_get_max_threads() > cl_thread_limit) {
#if _OPENMP
    int team_size = cl_thread_limit;
#else
    int team_size = 1;
#endif
    omp_set_num_threads(19);
    verify(__FILE__, __LINE__, 1);
#pragma omp parallel
    {
      verify(__FILE__, __LINE__, team_size);
      verify(__FILE__, __LINE__, team_size);
    }
    verify(__FILE__, __LINE__, 1);

    omp_set_nested(1);
#pragma omp parallel num_threads(3)
    {
      verify(__FILE__, __LINE__, 3);
#pragma omp master
#pragma omp parallel num_threads(21)
      {
        verify(__FILE__, __LINE__, team_size-2);
        verify(__FILE__, __LINE__, team_size-2);
      }
    }
    verify(__FILE__, __LINE__, 1);

    return failed;
  } else {
    fprintf(stderr, "This test is not applicable for max num_threads='%d'\n",
            omp_get_max_threads());
    return 0;

*/

void 
gms::common
::verify_thread_limit(const char * file_name,
                      const std::int32_t line_number,
                      const std::int32_t team_size)
{
  std::int32_t num_threads = omp_get_num_threads();
  if (team_size != num_threads) 
  {
#pragma omp critical(A)
    {
      char label[256];
      snprintf(label, sizeof(label), "%s:%d", file_name, line_number);
      failed = 1;
      std::printf("failed: %s: team_size(%d) != omp_get_num_threads(%d)\n",
                 label, team_size, num_threads);
    }
  }
}

#define XSTR(x) #x
#define STR(x) XSTR(x)

#define streqls(s1, s2) (!strcmp(s1, s2))

#define check(condition)                                                       \
  if (!(condition)) {                                                          \
    fprintf(stderr, "error: %s: %d: " STR(condition) "\n", __FILE__,           \
            __LINE__);                                                         \
    exit(1);                                                                   \
  }

__ATTR_ALWAYS_INLINE__ 
static inline 
std::int32_t get_integer()
{
   constexpr std::int32_t BUFFER_SIZE = 256;
   char buf[BUFFER_SIZE];
   size_t needed = omp_capture_affinity(buf, BUFFER_SIZE, NULL);
   check(needed < BUFFER_SIZE);
   n = sscanf(buf, "%d", &retval);
   check(n == 1);
   return retval;
}

__ATTR_ALWAYS_INLINE__
static inline 
char * get_string()
{
    constexpr std::int32_t BUFFER_SIZE = 256;
    char buf[BUFFER_SIZE];
    size_t needed = omp_capture_affinity(buf, BUFFER_SIZE, NULL);
    check(needed < BUFFER_SIZE);
    return strdup(buf);
}

__ATTR_ALWAYS_INLINE__
static inline
void check_integer(const char * formats[2], std::int32_t(*func)(),
                   const std::int32_t n_threads_outer,
                   const std::int32_t n_threads_inner)
{
    std::int32_t i;
    for (i = 0; i < 2; ++i) 
    {
    omp_set_affinity_format(formats[i]);
    #pragma omp parallel num_threads(n_threads_outer) // was 8
    {
      check(get_integer() == func());
      #pragma omp parallel num_threads(n_threads_inner) // was 3
      {
        check(get_integer() == func());
      }
      check(get_integer() == func());
    }
  }
}

void 
gms::common
::check_nesting_level(const std::int32_t n_threads_outer,
                      const std::int32_t n_threads_inner)
{
    const char* formats[2] = {"%{nesting_level}", "%L"};
    check_integer(formats, omp_get_level,n_threads_outer,n_threads_inner);
                   
}

void 
gms::common 
::check_thread_num(const std::int32_t n_threads_outer,
                   const std::int32_t n_threads_inner)
{
    const char* formats[2] = {"%{thread_num}", "%n"};
    check_integer(formats, omp_get_thread_num,n_threads_outer,n_threads_inner);
}

void 
gms::common
::check_num_threads(const std::int32_t n_threads_outer,
                    const std::int32_t n_threads_inner)
{
    const char* formats[2] = {"%{num_threads}", "%N"};
    check_integer(formats, omp_get_num_threads,n_threads_outer,n_threads_inner);
}

__ATTR_ALWAYS_INLINE__
static inline 
std::int32_t 
ancestor_helper() {return omp_get_ancestor_thread_num(omp_get_level() - 1); }

void 
gms::common 
::check_thread_ancestor(const std::int32_t n_threads_outer,
                        const std::int32_t n_threads_inner)
{
    const char* formats[2] = {"%{ancestor_tnum}", "%a"};
    check_integer(formats, ancestor_helper,n_threads_outer,n_threads_inner);
}

__ATTR_ALWAYS_INLINE__
static inline 
std::int32_t get_pid_helper() { return (static_cast<std::int32_t>(getpid()));}

void 
gms::common
::check_process_id(const std::int32_t n_threads_outer,
                   const std::int32_t n_threads_inner)
{
    const char* formats[2] = {"%{process_id}", "%P"};
    check_integer(formats, get_pid_helper,n_threads_outer,n_threads_inner);
}

void 
gms::common 
::check_host(const std::int32_t n_threads)
{
    typedef std::int32_t buffer_size_t;
    std::int32_t i;
    buffer_size_t buffer_size = 256;
    const char* formats[2] = {"%{host}", "%H"};
    char hostname[256];
    gethostname(hostname, buffer_size);
    for (i = 0; i < 2; ++i)
    {
        omp_set_affinity_format(formats[i]);
        #pragma omp parallel num_threads(n_threads)
        {
           char* host = get_string();
           check(streqls(host, hostname));
           std::free(host);
        }
    }
}
