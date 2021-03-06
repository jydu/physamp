//
// File: PhyloSampler.cpp
// Created by: Julien Dutheil
// Created on: Sunday, December 2nd 2007 16:48
//

/*
    This file is part of the PhySamp package.

    PhySamp is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    PhySamp is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
*/

// From the STL:
#include <iostream>
#include <iomanip>

using namespace std;

#include <Bpp/App/BppApplication.h>
#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Io/FileTools.h>
#include <Bpp/Text/TextTools.h>
#include <Bpp/Numeric/DataTable.h>
#include <Bpp/Numeric/Random/RandomTools.h>

// From bpp-seq:
#include <Bpp/Seq/SiteTools.h>
#include <Bpp/Seq/SequenceTools.h>
#include <Bpp/Seq/App/SequenceApplicationTools.h>

// From bpp-phyl:
#include <Bpp/Phyl/Tree.h>
#include <Bpp/Phyl/App/PhylogeneticsApplicationTools.h>
#include <Bpp/Phyl/Io/PhylipDistanceMatrixFormat.h>

using namespace bpp;

void help()
{
  (*ApplicationTools::message << "__________________________________________________________________________").endLine();
  (*ApplicationTools::message << "bppphysamp parameter1_name=parameter1_value").endLine();
  (*ApplicationTools::message << "      parameter2_name=parameter2_value ... param=option_file").endLine();
  (*ApplicationTools::message).endLine();
  (*ApplicationTools::message << "  Refer to the PhySamp manual for a list of available options.").endLine();
  (*ApplicationTools::message << "__________________________________________________________________________").endLine();
}

class Index {
  public:
    double distance;
    unsigned int i1, i2;

  public:
    Index(double dist, unsigned int i, unsigned int j) : distance(dist), i1(i), i2(j) {}

  public:
    bool operator==(const Index& index) const { return distance == index.distance; }
    bool operator<(const Index& index) const { return distance < index.distance; }
};

class Test {
  private:
    unsigned int pos_;

  public:
    Test(unsigned int pos) : pos_(pos) {}
    
  public:
    bool operator()(const Index& index) { return index.i1 == pos_ || index.i2 == pos_; }
};

int main(int args, char ** argv)
{
  cout << "******************************************************************" << endl;
  cout << "*           Bio++ Phylogenetic Sampler, version 1.1.0.           *" << endl;
  cout << "* Author: J. Dutheil                        Last Modif. 14/03/18 *" << endl;
  cout << "******************************************************************" << endl;
  cout << endl;
  
  if(args == 1)
  {
    help();
    return 0;
  }
  
  try {

  BppApplication bppphysamp(args, argv, "BppPhySamp");
  bppphysamp.startTimer();

  //Get sequences:
  Alphabet* alphabet      = SequenceApplicationTools::getAlphabet(bppphysamp.getParams());
  SequenceContainer* seqs = SequenceApplicationTools::getSequenceContainer(alphabet, bppphysamp.getParams());

  string inputMethod = ApplicationTools::getStringParameter("input.method", bppphysamp.getParams(), "tree");
  ApplicationTools::displayResult("Input method", inputMethod);

  unique_ptr< DistanceMatrix > dist;
  unique_ptr< TreeTemplate<Node> > tree;
  if(inputMethod == "tree")
  {
    tree.reset(dynamic_cast<TreeTemplate<Node> *>(PhylogeneticsApplicationTools::getTree(bppphysamp.getParams())));
    dist.reset(TreeTemplateTools::getDistanceMatrix(*tree));
    //PhylipDistanceMatrixFormat matIO;
    //matIO.write(*dist.get(), "matrix.txt");
  }
  else if(inputMethod == "matrix")
  {
    string distPath = ApplicationTools::getAFilePath("input.matrix", bppphysamp.getParams(), true, true);
    PhylipDistanceMatrixFormat matIO;
    dist.reset(matIO.readDistanceMatrix(distPath));
  }
  else throw Exception("Unknown input method: " + inputMethod);

  string deleteMeth = ApplicationTools::getStringParameter("deletion_method", bppphysamp.getParams(), "threshold");
  ApplicationTools::displayResult("Deletion method", deleteMeth);

  string critMeth = ApplicationTools::getStringParameter("choice_criterion", bppphysamp.getParams(), "length");
  ApplicationTools::displayResult("Sequence choice criterion", critMeth);

  //Compute lengths:
  vector<string> seqNames;
  vector<size_t> seqLen(dist->size());
  string name;
  for(size_t i = 0; i < dist->size(); i++)
  {
    name = dist->getName(i);
    if (critMeth == "length.complete")
      seqLen[i] = SequenceTools::getNumberOfCompleteSites(seqs->getSequence(name));
    else
      seqLen[i] = SequenceTools::getNumberOfSites(seqs->getSequence(name));
    seqNames.push_back(name);
  }

  //Sort matrix entries:
  vector<Index> distances;
  for (unsigned int i = 0; i < dist->size()-1; i++)
    for (unsigned int j = i+1; j < dist->size(); j++)
      distances.push_back(Index((*dist)(i, j), i , j));
  sort(distances.begin(), distances.end());

  if (deleteMeth == "random")
  {
    unsigned int sampleSize = ApplicationTools::getParameter<unsigned int>("sample_size", bppphysamp.getParams(), 10);
    ApplicationTools::displayResult("Sample size", sampleSize);
    vector<string> sample(sampleSize);
    RandomTools::getSample(seqNames, sample, false);
    seqNames = sample;
    
    double mini = -log(0.);
    for (unsigned int i =  0; i < seqNames.size() - 1; ++i)
      for (unsigned int j = i + 1; j < seqNames.size(); ++j)
      {
        double d = (*dist)(seqNames[i], seqNames[j]);
        if (d < mini) mini = d;
      }
    ApplicationTools::displayResult("Minimal distance in final data set:", mini);
  }
  else if (deleteMeth == "threshold")
  {
    double threshold = ApplicationTools::getDoubleParameter("threshold", bppphysamp.getParams(), 0.01);
    ApplicationTools::displayResult("Distance threshold", threshold);

    unsigned int rm = 0;
    while (distances[0].distance <= threshold)
    {
      //We need to chose between the two sequences:
      if (critMeth == "length" || critMeth == "length.complete")
      {
        if (seqLen[distances[0].i1] > seqLen[distances[0].i2]) rm = distances[0].i2;
        else rm = distances[0].i1;
      }
      else if (critMeth == "random")
      {
        if (RandomTools::flipCoin()) rm = distances[0].i2;
        else rm = distances[0].i1;
      }
      else throw Exception("Unknown criterion: " + critMeth);

      //Remove sequence in list:
      size_t pos = VectorTools::which(seqNames, dist->getName(rm));
      ApplicationTools::displayResult("Remove sequence", seqNames[pos]);
      seqNames.erase(seqNames.begin() + static_cast<ptrdiff_t>(pos)); 
        
      //Ignore all distances from this sequence:
      remove_if(distances.begin(), distances.end(), Test(rm));
      if (distances.size() == 0)
        throw Exception("Error, all sequences have been removed with this criterion!");
    }
    ApplicationTools::displayResult("Number of sequences kept:", seqNames.size());
  }
  else if (deleteMeth == "sample")
  {
    unsigned int sampleSize = ApplicationTools::getParameter<unsigned int>("sample_size", bppphysamp.getParams(), 10);
    ApplicationTools::displayResult("Sample size", sampleSize);
    
    unsigned int rm = 0;
    while (seqNames.size() > sampleSize)
    {
      //We need to chose between the two sequences:
      if (critMeth == "length" || critMeth == "length.complete")
      {
        if (seqLen[distances[0].i1] > seqLen[distances[0].i2]) rm = distances[0].i2;
        else rm = distances[0].i1;
      }
      else if (critMeth == "random")
      {
        if (RandomTools::flipCoin()) rm = distances[0].i2;
        else rm = distances[0].i1;
      }
      else throw Exception("Unknown criterion: " + critMeth);

      //Remove sequence in list:
      size_t pos = VectorTools::which(seqNames, dist->getName(rm));
      ApplicationTools::displayResult("Remove sequence", seqNames[pos]);
      seqNames.erase(seqNames.begin() + static_cast<ptrdiff_t>(pos)); 
        
      //Ignore all distances from this sequence:
      remove_if(distances.begin(), distances.end(), Test(rm));
    }
    ApplicationTools::displayResult("Minimal distance in final data set:", distances[0].distance);
  }
  else throw Exception("Unknown deletion method: " + deleteMeth + ".");

  //Write sequences to file:
  AlignedSequenceContainer asc(alphabet);
  for (size_t i = 0; i < seqNames.size(); i++)
    asc.addSequence(seqs->getSequence(seqNames[i]));
   
  SequenceApplicationTools::writeAlignmentFile(asc, bppphysamp.getParams());

  //Write tree file:
  if (ApplicationTools::getStringParameter("output.tree.file", bppphysamp.getParams(), "None") != "None") {
    vector<string> allSeqNames(seqs->getSequencesNames());
    vector<string> removedSeqNames;
    VectorTools::diff(allSeqNames, seqNames, removedSeqNames);
    for (size_t i = 0; i < removedSeqNames.size(); ++i) {
      TreeTemplateTools::dropLeaf(*tree, removedSeqNames[i]);
    }
    PhylogeneticsApplicationTools::writeTree(*tree, bppphysamp.getParams(), "output.", "", true, true, false);
  }

  bppphysamp.done();
  }
  catch (exception& e)
  {
    cout << endl;
    cout << "_____________________________________________________" << endl;
    cout << "ERROR!!!" << endl;
    cout << e.what() << endl;
    return 1;
  }

  return 0;
}

