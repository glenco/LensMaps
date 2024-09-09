#include <slsimlib.h>
#include <gridmap.h>
#include <ctime>

int main(int argc, char** argv)
{
  std::time_t t0;
  std::time(&t0);
  
  long seed = -1234567890;
  int gridsize = 1024;
  std::string outputfile = "example_output";
  
  const COSMOLOGY cosmo(CosmoParamSet::Planck18);
  Lens lens(&seed,10,cosmo);
  
  std::cout
  << std::endl;
  
  std::cout
  << "cosmology" << std::endl
  << "  h: " << cosmo.gethubble() << std::endl
  << "  Omega_m: " << cosmo.getOmega_matter() << std::endl
  << "  Omega_L: " << cosmo.getOmega_lambda() << std::endl;
  
  std::vector<double> source_redshifts = {0.5,1.0,1.5,2.0,2.5};
  
  // read lens plane file names
  std::vector<std::string> mass_files;
  {
    std::ifstream file("example_planes.txt");
    if (!file.is_open()){
      std::cerr << "file example_planes.txt cann't be opened." << std::endl;
      throw std::runtime_error("no file");
    }
    std::string line;
    while(std::getline(file,line)){
      if(line[0] != '#') mass_files.push_back(line);
    }
  }
  
  Point_2d center(0,0);
  double range = 0;
  // read mass planes and insert them into the lens
  for(std::string filename : mass_files){
    LensHaloMassMap map(filename
                  ,PixelMapType::pix_map
                  ,0                        // pad sides ??
                  ,true                     // subtract mean
                  ,cosmo
                  );
    
    // this is for if the redshif and mass conversion are not given in the header of the fits file
    //LensHaloMassMap map(filename
    //              ,PixelMapType::pix_map
    //              ,massconvertion    // convertion factor from pixel units to solar masses
    //              ,redshift          // redshift of mass plane
    //              ,0                 // pad sides ??
    //              ,true              // subtract mean
    //              ,cosmo
    //              );
    

    
    std::cout
    << "lens map " << filename << std::endl
    << "  size: " << map.getNx() << "x" << map.getNy() << std::endl
    << "  redshift: " << map.getZlens() << std::endl
    << "  range : " << map.getRangeMpc() << " Mpc = " << map.getRangeRad()/arcminTOradians << " arcmin" << std::endl;
    
    center += map.getCenter()/cosmo.angDist(map.getZlens());  // convert to angle
    range = MAX(map.getRangeRad(),range);
    lens.moveinMainHalo(map,true);
  }
  center /= mass_files.size();
  
  std::cout
  << "lens" << std::endl
  << "  number of planes: " << lens.getNplanes() << std::endl
  << "  source redshifts: ";
  for(auto z: source_redshifts)
    std::cout << z << " ";
  std::cout << std::endl;

  std::cout
  << "grid" << std::endl
  << "  size: " << gridsize << std::endl
  << "  center: (" << center[0] << ", " << center[1] << ")" << std::endl
  << "  range: " << (range/degreesTOradians) << " deg" << std::endl;
  
  // loop through redshifts
  for(auto z: source_redshifts)
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
    GridMap grid(&lens, gridsize, center.x , range);
    std::cout << "done" << std::endl;
    
    // write kappa map
    std::cout << "  kappa " << std::flush;
    grid.writeFits<float>(LensingVariable::KAPPA,"!" + outfile);
    //grid.writeFitsUniform(mmap->getCenter(), mapsize, mapsize, KAPPA, "!" + outfile);
    std::cout << "done" << std::endl;
    
    // write gamma map
    std::cout << "  gamma " << std::flush;
    grid.writeFits<float>(LensingVariable::GAMMA1,"!" + outfile);
    grid.writeFits<float>(LensingVariable::GAMMA2,"!" + outfile);
    //grid.writeFitsUniform(mmap->getCenter(), mapsize, mapsize, GAMMA1, "!" + outfile);
    //grid.writeFitsUniform(mmap->getCenter(), mapsize, mapsize, GAMMA2, "!" + outfile);
    std::cout << "done" << std::endl;
  }
  
  std::time_t t1;
  std::time(&t1);
  
  std::cout
  << std::endl
  << "finished in " << (std::difftime(t1,t0)/60.) << " mins" << std::endl;
  return EXIT_SUCCESS;
}
