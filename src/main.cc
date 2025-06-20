// This file is part of monofonIC (MUSIC2)
// A software package to generate ICs for cosmological simulations
// Copyright (C) 2020 by Oliver Hahn
// 
// monofonIC is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// monofonIC is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include <cmath>
#include <complex>
#include <iostream>
#include <fstream>
#include <thread>
#include <cfenv>

#if defined(_OPENMP)
#include <omp.h>
#endif

#include <general.hh>
#include <ic_generator.hh>
#include <cosmology_parameters.hh>
#include <particle_plt.hh>


// initialise with "default" values
namespace CONFIG{
int  MPI_thread_support = -1;
int  MPI_task_rank = 0;
int  MPI_task_size = 1;
bool MPI_ok = false;
bool MPI_threads_ok = false;
bool FFTW_threads_ok = false;
int  num_threads = 1;
}

size_t global_mem_high_mark, local_mem_high_mark;

#include "system_stat.hh"
#include "memory_stat.hh"


#include <exception>
#include <stdexcept>
 
void handle_eptr(std::exception_ptr eptr) // passing by value is ok
{
    try {
        if (eptr) {
            std::rethrow_exception(eptr);
        }
    } catch(const std::exception& e) {
        music::elog << "This happened: \"" << e.what() << "\"" << std::endl;
    }
}

/**
 * @brief the main routine of MUSIC2-monofonIC
 * 
 * @param argc 
 * @param argv 
 * @return int 
 */
int main( int argc, char** argv )
{

#if defined(NDEBUG)
    music::logger::set_level(music::log_level::info);
#else
    music::logger::set_level(music::log_level::debug);
#endif

    global_mem_high_mark = local_mem_high_mark = 0;

    //------------------------------------------------------------------------------
    // initialise MPI 
    //------------------------------------------------------------------------------
    
#if defined(USE_MPI)
    int thread_wanted = MPI_THREAD_MULTIPLE; // MPI_THREAD_FUNNELED
    MPI_Init_thread(&argc, &argv, thread_wanted, &CONFIG::MPI_thread_support);
    CONFIG::MPI_threads_ok = CONFIG::MPI_thread_support >= thread_wanted;
    MPI_Comm_rank(MPI_COMM_WORLD, &CONFIG::MPI_task_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &CONFIG::MPI_task_size);
    CONFIG::MPI_ok = true;

    // set up lower logging levels for other tasks
    if( CONFIG::MPI_task_rank!=0 )
    {
        music::logger::set_level(music::log_level::error);
    }
#endif

    // Ascii ART logo. generated via http://patorjk.com/software/taag/#p=display&f=Nancyj&t=monofonIC
    music::ilog << "\n"
                << " The unigrid version of MUSIC-2         .8888b                   dP  a88888b. \n"
                << "                                        88   \"                   88 d8\'   `88 \n"
                << "  88d8b.d8b. .d8888b. 88d888b. .d8888b. 88aaa  .d8888b. 88d888b. 88 88        \n"
                << "  88\'`88\'`88 88\'  `88 88\'  `88 88\'  `88 88     88\'  `88 88\'  `88 88 88        \n"
                << "  88  88  88 88.  .88 88    88 88.  .88 88     88.  .88 88    88 88 Y8.   .88 \n"
                << "  dP  dP  dP `88888P\' dP    dP `88888P\' dP     `88888P\' dP    dP dP  Y88888P\' \n" << std::endl;

    // git and versioning info:
    music::ilog << "Version: git rev.: " << GIT_REV << ", tag: " << GIT_TAG << ", branch: " << GIT_BRANCH << std::endl;
    
    // Compilation CMake configuration, time etc info:
    music::ilog << "This " << CMAKE_BUILDTYPE_STR << " build was compiled at " << __TIME__ << " on " <<  __DATE__ << std::endl;

#ifdef __GNUC__
    music::ilog << "Compiled with GNU C++ version " << __VERSION__ <<std::endl;
#else
    music::ilog << "Compiled with " << __VERSION__ << std::endl;
#endif

    
    music::ilog << "-------------------------------------------------------------------------------" << std::endl;
    music::ilog << "Compile time options : " << std::endl;
    music::ilog << "                       Precision : " << CMAKE_PRECISION_STR << std::endl;
    music::ilog << "               Convolutions (LPT): " << CMAKE_CONVOLVER_STR << std::endl;
    music::ilog << "               Convolutions (PNG): " << CMAKE_CONVOLVER_PNG_STR << std::endl;
    music::ilog << "                             PLT : " << CMAKE_PLT_STR << std::endl;
    music::ilog << "-------------------------------------------------------------------------------" << std::endl;


    //------------------------------------------------------------------------------
    // Parse command line options
    //------------------------------------------------------------------------------

    if (argc != 2)
    {
        // print_region_generator_plugins();
        cosmology::print_ParameterSets();
        print_TransferFunction_plugins();
        print_RNG_plugins();
        print_output_plugins();

        music::elog << "In order to run, you need to specify a parameter file!\n" << std::endl;
        exit(0);
    }

    // open the configuration file 
    config_file the_config(argv[1]);
    std::string log_filename = the_config.get_path_relative_to_config("log.txt");
    music::logger::set_output(log_filename);
    music::ilog << "                         argv[1] : " << argv[1] << std::endl;
    music::ilog << "                 config_basename : " << the_config.get_path_relative_to_config("") << std::endl;
    music::ilog << "                        log file : " << log_filename << std::endl;
    music::ilog << "-------------------------------------------------------------------------------" << std::endl;

    //------------------------------------------------------------------------------
    // Set up FFTW
    //------------------------------------------------------------------------------

#if defined(USE_FFTW_THREADS)
  #if defined(USE_MPI)
    if (CONFIG::MPI_threads_ok)
        CONFIG::FFTW_threads_ok = FFTW_API(init_threads)();
  #else
    CONFIG::FFTW_threads_ok = FFTW_API(init_threads)();
  #endif 
#endif

#if defined(USE_MPI)
    FFTW_API(mpi_init)();
#endif

    CONFIG::num_threads = the_config.get_value_safe<unsigned>("execution", "NumThreads",std::thread::hardware_concurrency());
    
#if defined(USE_FFTW_THREADS)
    if (CONFIG::FFTW_threads_ok)
        FFTW_API(plan_with_nthreads)(CONFIG::num_threads);
#endif

    //------------------------------------------------------------------------------
    // Set up OpenMP
    //------------------------------------------------------------------------------

#if defined(_OPENMP)
    omp_set_num_threads(CONFIG::num_threads);
#endif

    std::feclearexcept(FE_ALL_EXCEPT);

    //------------------------------------------------------------------------------
    // Write code configuration to screen
    //------------------------------------------------------------------------------
    // hardware related infos
    music::ilog << std::setw(32) << std::left << "CPU vendor string" << " : " << SystemStat::Cpu().get_CPUstring() << std::endl;
    
    // multi-threading related infos
    music::ilog << std::setw(32) << std::left << "Available HW threads / task" << " : " << std::thread::hardware_concurrency() << " (" << CONFIG::num_threads << " used)" << std::endl;

    // memory related infos
    SystemStat::Memory mem;

    unsigned availpmem = mem.get_AvailMem()/1024/1024;
    unsigned usedpmem = mem.get_UsedMem()/1024/1024;
    unsigned maxpmem = availpmem, minpmem = availpmem;
    unsigned maxupmem = usedpmem, minupmem = usedpmem;
    
#if defined(USE_MPI)
    unsigned temp = 0;
    MPI_Allreduce(&minpmem,&temp,1,MPI_UNSIGNED,MPI_MIN,MPI_COMM_WORLD);  minpmem = temp;
    MPI_Allreduce(&maxpmem,&temp,1,MPI_UNSIGNED,MPI_MAX,MPI_COMM_WORLD);  maxpmem = temp;
    MPI_Allreduce(&minupmem,&temp,1,MPI_UNSIGNED,MPI_MIN,MPI_COMM_WORLD); minupmem = temp;
    MPI_Allreduce(&maxupmem,&temp,1,MPI_UNSIGNED,MPI_MAX,MPI_COMM_WORLD); maxupmem = temp;
#endif
    music::ilog << std::setw(32) << std::left << "Total system memory (phys)" << " : " << mem.get_TotalMem()/1024/1024 << " Mb" << std::endl;
    music::ilog << std::setw(32) << std::left << "Used system memory (phys)" << " : " << "Max: " << maxupmem << " Mb, Min: " << minupmem << " Mb" << std::endl;
    music::ilog << std::setw(32) << std::left << "Available system memory (phys)" << " : " <<  "Max: " << maxpmem << " Mb, Min: " << minpmem << " Mb" << std::endl;
    
    // MPI related infos
#if defined(USE_MPI)
    music::ilog << std::setw(32) << std::left << "MPI is enabled" << " : " << "yes (" << CONFIG::MPI_task_size << " tasks)" << std::endl;
    music::dlog << std::setw(32) << std::left << "MPI version" << " : " << MPI::get_version() << std::endl;
#else
    music::ilog << std::setw(32) << std::left << "MPI is enabled" << " : " << "no" << std::endl;
#endif
    music::ilog << std::setw(32) << std::left << "MPI supports multi-threading" << " : " << (CONFIG::MPI_threads_ok? "yes" : "no") << std::endl;
    
    // Kernel related infos
    SystemStat::Kernel kern;
    auto kinfo = kern.get_kernel_info();
    music::ilog << std::setw(32) << std::left << "OS/Kernel version" << " : " << kinfo.kernel << " version " << kinfo.major << "." << kinfo.minor << " build " << kinfo.build_number << std::endl;

    // FFTW related infos
    music::ilog << std::setw(32) << std::left << "FFTW version" << " : " << FFTW_API(version) << std::endl;
    music::ilog << std::setw(32) << std::left << "FFTW supports multi-threading" << " : " << (CONFIG::FFTW_threads_ok? "yes" : "no") << std::endl;
    music::ilog << std::setw(32) << std::left << "FFTW mode" << " : ";
#if defined(FFTW_MODE_PATIENT)
	music::ilog << "FFTW_PATIENT" << std::endl;
#elif defined(FFTW_MODE_MEASURE)
    music::ilog << "FFTW_MEASURE" << std::endl;
#else
	music::ilog << "FFTW_ESTIMATE" << std::endl;
#endif

    ///////////////////////////////////////////////////////////////////////
    // Initialise plug-ins
    try
    {
        ic_generator::initialise( the_config );
    }catch(...){
        handle_eptr( std::current_exception() );
        music::elog << "Problem during initialisation. See error(s) above. Exiting..." << std::endl;
        #if defined(USE_MPI) 
        MPI_Finalize();
        #endif
        return 1;
    }
    ///////////////////////////////////////////////////////////////////////


    ///////////////////////////////////////////////////////////////////////
    // do the job...
    ic_generator::run( the_config );
    ///////////////////////////////////////////////////////////////////////


    ///////////////////////////////////////////////////////////////////////
    // call the destructor of plugins before tearing down MPI
    ic_generator::reset();
    ///////////////////////////////////////////////////////////////////////

    music::ilog << "-------------------------------------------------------------------------------" << std::endl;
    size_t peak_mem = memory::getPeakRSS();
#if defined(USE_MPI)
    size_t peak_mem_max{0};
    MPI_Allreduce(&peak_mem, &peak_mem_max, 1, MPI_UNSIGNED_LONG_LONG, MPI_MAX, MPI_COMM_WORLD);
    peak_mem = peak_mem_max;
#endif

    if( peak_mem > (1ull<<30) )
        music::ilog << "Peak memory usage was " << peak_mem /(1ull<<30) << " GBytes / task" << std::endl;
    else 
        music::ilog << "Peak memory usage was " << peak_mem /(1ull<<20) << " MBytes / task" << std::endl;


#if defined(USE_MPI)
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
#endif
        
    music::ilog << "Done. Have a nice day!\n" << std::endl;

    return 0;
}
