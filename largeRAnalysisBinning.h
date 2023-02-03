#ifndef MY_HEADER_FILE_H
#define MY_HEADER_FILE_H

#include <vector>
  namespace {
    //R=0.4
    //Truth Binning
    std::vector<std::vector<double>> binsR4_Truth;
  


  //R=0.4
  //Truth Binning
  //vector<std::vector<double>> binsR4_Truth;
  binsR4_Truth.reserve(8);
  binsR4_Truth.emplace_back(std::vector<double>{66,83,100, 117, 138, 158, 178, 199, 220, 241, 262, 294, 336, 384, 439, 502, 573, 656, 750, 857, 980, 1120,1300}); //0-10%
  binsR4_Truth.emplace_back(std::vector<double>{66,83,100, 117, 138, 158, 178, 199, 220, 241, 262, 294, 336, 384, 439, 502, 573, 656, 750, 857, 980, 1120, 1300}); //10-20%
  binsR4_Truth.emplace_back(std::vector<double>{75,92,109, 126, 143, 160, 177, 194, 211, 228, 245, 262, 294, 336, 384, 439, 502, 573, 656, 750, 857, 980, 1120}); //20-30%
  binsR4_Truth.emplace_back(std::vector<double>{66,83,100, 117, 138, 158, 177, 198, 230, 262, 294, 336, 384, 439, 502, 573, 656, 750, 857, 1000});//30-40%
  binsR4_Truth.emplace_back(std::vector<double>{66,83,100, 117, 138, 158, 177, 198, 230, 262, 294, 336, 384, 439, 502, 573, 656, 750, 857, 1000}); //40-50%
  binsR4_Truth.emplace_back(std::vector<double>{66,83,100, 117, 138, 158, 177, 198, 220, 241, 262, 294, 336, 384, 439, 502, 573, 656, 750, 857, 1000}); //50-60%
  binsR4_Truth.emplace_back(std::vector<double>{66,83,100, 117, 138, 158, 177, 197, 220, 241, 262, 294, 336, 384, 439, 502, 573, 656, 750, 1000}); //60-70%
  binsR4_Truth.emplace_back(std::vector<double>{66,83,100, 117, 138, 158, 177, 197, 220, 241, 262, 294, 336, 384, 439, 502, 573, 656, 750, 1000}); //70-80%
  //Reco Binning
  std::vector<std::vector<double>> binsR4;
  binsR4.reserve(8);
  binsR4.emplace_back(std::vector<double>{100, 117, 138, 158, 178, 199, 220, 241, 262, 294, 336, 384, 439, 502, 573, 656, 750, 857, 980, 1120,1300}); //0-10%
  binsR4.emplace_back(std::vector<double>{100, 117, 138, 158, 178, 199, 220, 241, 262, 294, 336, 384, 439, 502, 573, 656, 750, 857, 980, 1120, 1300}); //10-20%
  binsR4.emplace_back(std::vector<double>{109, 126, 143, 160, 177, 194, 211, 228, 245, 262, 294, 336, 384, 439, 502, 573, 656, 750, 857, 980, 1120}); //20-30%
  binsR4.emplace_back(std::vector<double>{100, 117, 138, 158, 177, 198, 230, 262, 294, 336, 384, 439, 502, 573, 656, 750, 857, 1000}); //30-40%
  binsR4.emplace_back(std::vector<double>{100, 117, 138, 158, 177, 198, 230, 262, 294, 336, 384, 439, 502, 573, 656, 750, 857, 1000}); //40-50%
  binsR4.emplace_back(std::vector<double>{100, 117, 138, 158, 177, 198, 220, 241, 262, 294, 336, 384, 439, 502, 573, 656, 750, 857, 1000}); //50-60%
  binsR4.emplace_back(std::vector<double>{100, 117, 138, 158, 177, 197, 220, 241, 262, 294, 336, 384, 439, 502, 573, 656, 750, 1000}); //60-70%
  binsR4.emplace_back(std::vector<double>{100, 117, 138, 158, 177, 197, 220, 241, 262, 294, 336, 384, 439, 502, 573, 656, 750, 1000}); //70-80~

  std::vector<std::vector<double>> binsR10_Truth;
  binsR10_Truth.reserve(8);
  binsR10_Truth.emplace_back(std::vector<double>{241, 262, 294, 336, 384, 439, 502, 573, 656, 750, 857, 980, 1120,1300}); //0-10%
  binsR10_Truth.emplace_back(std::vector<double>{220 ,241, 262, 294, 336, 384, 439, 502, 573, 656, 750, 857, 980, 1120, 1300}); //10-20%
  binsR10_Truth.emplace_back(std::vector<double>{211, 228, 245, 262, 294, 336, 384, 439, 502, 573, 656, 750, 857, 980, 1120}); //20-30%
  binsR10_Truth.emplace_back(std::vector<double>{177, 198, 230, 262, 294, 336, 384, 439, 502, 573, 656, 750, 857, 1000}); //30-40%
  binsR10_Truth.emplace_back(std::vector<double>{177,198, 230, 262, 294, 336, 384, 439, 502, 573, 656, 750, 857, 1000}); //40-50%
  binsR10_Truth.emplace_back(std::vector<double>{177, 198, 220, 241, 262, 294, 336, 384, 439, 502, 573, 656, 750, 857, 1000}); //50-60%
  binsR10_Truth.emplace_back(std::vector<double>{177, 197, 220, 241, 262, 294, 336, 384, 439, 502, 573, 656, 750, 1000}); //60-70%
  binsR10_Truth.emplace_back(std::vector<double>{177, 197, 220, 241, 262, 294, 336, 384, 439, 502, 573, 656, 750, 1000}); //70-80%



  //Reco Binning
  std::vector<std::vector<double>> binsR10;
  binsR10.reserve(8);
  binsR10.emplace_back(std::vector<double>{294, 336, 384, 439, 502, 573, 656, 750, 857, 980, 1120,1300}); //0-10%
  binsR10.emplace_back(std::vector<double>{262, 294, 336, 384, 439, 502, 573, 656, 750, 857, 980, 1120, 1300}); //10-20%
  binsR10.emplace_back(std::vector<double>{245, 262, 294, 336, 384, 439, 502, 573, 656, 750, 857, 980, 1120}); //20-30%
  binsR10.emplace_back(std::vector<double>{230, 262, 294, 336, 384, 439, 502, 573, 656, 750, 857, 1000}); //30-40%
  binsR10.emplace_back(std::vector<double>{230, 262, 294, 336, 384, 439, 502, 573, 656, 750, 857, 1000}); //40-50%
  binsR10.emplace_back(std::vector<double>{220, 241, 262, 294, 336, 384, 439, 502, 573, 656, 750, 857, 1000}); //50-60%
  binsR10.emplace_back(std::vector<double>{220, 241, 262, 294, 336, 384, 439, 502, 573, 656, 750, 1000}); //60-70%
  binsR10.emplace_back(std::vector<double>{220, 241, 262, 294, 336, 384, 439, 502, 573, 656, 750, 1000}); //70-80%
  }
#endif
