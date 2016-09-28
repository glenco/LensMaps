#include <slsimlib.h>
#include <gridmap.h>
#include <ctime>

int main(int argc, char** argv)
{
    std::time_t t0;
    std::time(&t0);
    
    std::string paramfile = (argc > 1) ? argv[1] : "paramfile";
    
    long seed = -1234567890;
    
    InputParams params = InputParams(paramfile);
    
    std::string outputfile;
    params.get("outputfile", outputfile);
    
    Lens lens(params, &seed);
    const COSMOLOGY& cosmo = lens.getCosmo();
    
    std::cout
    << std::endl;
     
    std::cout
    << "cosmology" << std::endl
    << "  h: " << cosmo.gethubble() << std::endl
    << "  Omega_m: " << cosmo.getOmega_matter() << std::endl
    << "  Omega_L: " << cosmo.getOmega_lambda() << std::endl;
    
    std::vector<double> redshifts;
    std::string redshifts_string;
    
    if(params.get("redshifts", redshifts_string))
    {
        std::stringstream sstr(redshifts_string);
        double z;
        char d;
        while(sstr.good())
        {
            sstr >> z;
            redshifts.push_back(z);
            sstr >> d;
        }
    }
    else
    {
        redshifts.push_back(lens.getSourceZ());
    }
    
    std::cout
    << "lens" << std::endl
    << "  number of planes: " << lens.getNplanes() << std::endl
    << "  source redshifts: ";
    for(auto z: redshifts)
        std::cout << z << " ";
    std::cout << std::endl;
    
    LensHaloMassMap* mmap = lens.getMainHalo<LensHaloMassMap>(0);
    
    std::cout
    << "mass map" << std::endl
    << "  size: " << mmap->getN() << std::endl;
    
    std::size_t gridsize = std::pow(2, 1 + std::ceil(log2(mmap->getN())));
    
    std::cout
    << "grid" << std::endl
    << "  size: " << gridsize << std::endl
    << "  center: (" << mmap->getCenter()[0] << ", " << mmap->getCenter()[1] << ")" << std::endl
    << "  range: " << (mmap->getRangeRad()/pi*180) << " deg" << std::endl;
    
    // we can ask for an higher resolution map 
    //int n = 512;
    int n = mmap->getN();
    
    std::cout
    << "lensing map" << std::endl
    << "  resolution: " << n << "x" << n << std::endl;
    
    // loop through redshifts
    for(auto z: redshifts)
    {
        std::stringstream sstr;
        sstr << outputfile;
        sstr << ".z";
        sstr << z;
        std::string outfile = sstr.str();
        
        std::cout << "source redshift " << z << std::endl;
        
        // set the source plane redshift
        lens.ResetSourcePlane(z, false);
        
        // create a fixed grid
        std::cout << "  grid " << std::flush; 
        GridMap grid(&lens, gridsize, mmap->getCenter(), mmap->getRangeRad(), mmap->getRangeRad());
        std::cout << "done" << std::endl;
        
        // write kappa map
        std::cout << "  kappa " << std::flush;
        grid.writeFitsUniform(mmap->getCenter(), n, n, KAPPA, "!" + outfile);
        std::cout << "done" << std::endl;
        
        // write gamma map
        std::cout << "  gamma " << std::flush;
        grid.writeFitsUniform(mmap->getCenter(), n, n, GAMMA1, "!" + outfile);
        grid.writeFitsUniform(mmap->getCenter(), n, n, GAMMA2, "!" + outfile);
        std::cout << "done" << std::endl;
    }
    
    std::time_t t1;
    std::time(&t1);
    
    std::cout
    << std::endl
    << "finished in " << (std::difftime(t1,t0)/60.) << " mins" << std::endl;
    
    return EXIT_SUCCESS;
}
