//
// File: AlignmentOptimizer.cpp
// Created by: Julien Dutheil
// Created on: October 24th 2014 11:15
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
#include <Bpp/Text/KeyvalTools.h>
#include <Bpp/Numeric/DataTable.h>
#include <Bpp/Numeric/Random/RandomTools.h>

// From bpp-seq:
#include <Bpp/Seq/Alphabet/AlphabetTools.h>
#include <Bpp/Seq/Container/SequenceContainerTools.h>
#include <Bpp/Seq/SiteTools.h>
#include <Bpp/Seq/App/SequenceApplicationTools.h>

// From bpp-phyl:
#include <Bpp/Phyl/TreeTemplate.h>
#include <Bpp/Phyl/Distance/HierarchicalClustering.h>
#include <Bpp/Phyl/App/PhylogeneticsApplicationTools.h>

using namespace bpp;

void help()
{
  (*ApplicationTools::message << "__________________________________________________________________________").endLine();
  (*ApplicationTools::message << "bppalnoptim parameter1_name=parameter1_value").endLine();
  (*ApplicationTools::message << "      parameter2_name=parameter2_value ... param=option_file").endLine();
  (*ApplicationTools::message).endLine();
  (*ApplicationTools::message << "  Refer to the PhySamp manual for a list of available options.").endLine();
  (*ApplicationTools::message << "__________________________________________________________________________").endLine();
}

//Note: this is too slow!
double overlapDistance(const Sequence& seq1, const Sequence& seq2) {
  //double n1 = 0, n2 = 0, n12 = 0;
  double n12 = 0;
  const Alphabet* a = seq1.getAlphabet();
  for (size_t i = 0; i < seq1.size(); ++i) {
    //if (!a->isGap(seq1[i])) n1++;
    //if (!a->isGap(seq2[i])) n2++;
    if (!a->isGap(seq1[i]) && !a->isGap(seq2[i])) n12++;
  }
  //return (1. - 2. * n12 / (n1 + n2)); 
  return (static_cast<double>(seq1.size()) - n12);
}

double overlapDistance(const vector<bool>& seq1, const vector<bool>& seq2) {
  double n12 = 0;
  for (size_t i = 0; i < seq1.size(); ++i) {
    if (seq1[i] && seq2[i]) n12++;
  }
  return (static_cast<double>(seq1.size()) - n12);
}

struct Group {
  public:
    int id;
    char direction;
    unsigned int nbSites;
    unsigned int nbSequences;
    unsigned int nbChars;

  public:
    Group(int nodeId,
        char dir,
        unsigned int nsi,
        unsigned int nse,
        unsigned int nch):
      id(nodeId), direction(dir),
      nbSites(nsi), nbSequences(nse), nbChars(nch)
    {}
};

class AlignmentPartitionScores {
  public:
    int id; //Needed later when in a vector and not a map
    unsigned int nbSitesUp;
    unsigned int nbSitesDown;
    unsigned int nbSequencesUp;
    unsigned int nbSequencesDown;
    unsigned int nbCharsUp;
    unsigned int nbCharsDown;

  public:
    AlignmentPartitionScores():
      id(-1),
      nbSitesUp(0),
      nbSitesDown(0),
      nbSequencesUp(0),
      nbSequencesDown(0),
      nbCharsUp(0),
      nbCharsDown(0)
    {}
    
    virtual ~AlignmentPartitionScores() {}

  public:
    Group getGroupUp() const {
      return Group(id, 'u', nbSitesUp, nbSequencesUp, nbCharsUp);
    }

    Group getGroupDown() const {
      return Group(id, 'd', nbSitesDown, nbSequencesDown, nbCharsDown);
    }
};

//Complete site version:
class AlignmentPartitionScores1:
  public AlignmentPartitionScores
{
  public:
    vector<bool> up;
    vector<bool> down;

  public:
    AlignmentPartitionScores1():
      up(), down() {}
    virtual ~AlignmentPartitionScores1() {}
};

//Version allowing a certain amount of gaps:
class AlignmentPartitionScores2:
  public AlignmentPartitionScores
{
  public:
    vector<unsigned int> up;
    vector<unsigned int> down;
  
  public:
    AlignmentPartitionScores2():
      up(), down() {}
    virtual ~AlignmentPartitionScores2() {}
};

class Comparator {
  private:
    char method_;

  public:
    Comparator(char method): method_(method) {}

  public:
    bool operator() (Group a, Group b) {
      if (method_ == 'a') {
        if (a.nbSequences == b.nbSequences) {
          if (a.nbSites == b.nbSites) {
            return (a.nbChars > b.nbChars);
          } else {
            return (a.nbSites > b.nbSites);
          }
        } else {
          return (a.nbSequences > b.nbSequences);
        }
      } else if (method_ == 'b') {
        return (a.nbChars > b.nbChars);
      } else {
        return (true);
      }
    }

    bool isImproving(Group& a, size_t nbSites, size_t nbChars) {
      if (method_ == 'a')
        return (a.nbSites > nbSites);
      else if (method_ == 'b') {
        return (a.nbChars > nbChars);
      } else {
        return false;
      }
    }
};


#define Index2 map<int, AlignmentPartitionScores2>

// Recursive function.
// This will initialize the scores map and set the down variables.
// For convenience, we return the scores for the root subtree.
AlignmentPartitionScores2& gapOptimizerFirstTreeTraversal2(const Node& node, const SiteContainer& sites, Index2& scores, double threshold, bool filterGaps, bool filterUnresolved, const vector<unsigned int>& reference, size_t nbSequencesRef)
{
  AlignmentPartitionScores2& nodeScores = scores[node.getId()]; //This eventually creates scores for this node.
  size_t nbSites = sites.getNumberOfSites();
  nodeScores.id = node.getId();
  nodeScores.down.assign(nbSites, 0);
  nodeScores.nbSitesDown = 0;
  nodeScores.nbCharsDown = 0;
  
  //We distinguish two cases, whether the node is a leaf or an inner node:
  if (node.isLeaf()) {
    //This is a leaf, we initialize arrays from the sequence data.
    const Sequence& seq = sites.getSequence(node.getName()); //We assume that the alignment contains all leaves of the tree.
    //This one at least is trivial:
    nodeScores.nbSequencesDown = 1;
    //Now we need to initialize the bit vector and count sites:
    for (size_t i = 0; i < nbSites; ++i) {
      bool test = !((filterGaps && seq.getAlphabet()->isGap(seq[i])) || (filterUnresolved && seq.getAlphabet()->isUnresolved(seq[i])));
      nodeScores.down[i] += (test ? 1 : 0);
    }
  } else {
    //This is an inner node.
    //All calculations are performed recursively from the results of son nodes:
    nodeScores.nbSequencesDown = 0;
    for (size_t k = 0; k < node.getNumberOfSons(); ++k) {
      //We loop over all son nodes:
      const Node* son = node.getSon(k);
      //We need to make computations for the subtree first and get the results:
      AlignmentPartitionScores2& sonScores = gapOptimizerFirstTreeTraversal2(*son, sites, scores, threshold, filterGaps, filterUnresolved, reference, nbSequencesRef);
      //Now we make computations for this node:
      nodeScores.nbSequencesDown += sonScores.nbSequencesDown;
      for (size_t i = 0; i < nbSites; ++i)
        nodeScores.down[i] += sonScores.down[i];
    }
  }
  //Finally we count all sites:
  for (size_t i = 0; i < nbSites; ++i) {
    if (static_cast<double>(nodeScores.down[i] + reference[i]) / static_cast<double>(nodeScores.nbSequencesDown + nbSequencesRef) >= threshold) {
      nodeScores.nbSitesDown++;
      nodeScores.nbCharsDown += nodeScores.down[i];
    }
  }
  return nodeScores;
}

// Recursive function.
// This will initialize the scores map and set the up variables.
// This must be ran after the first tree traversal, so that all scores are already initialized.
void gapOptimizerSecondTreeTraversal2(const Node& node, const SiteContainer& sites, Index2& scores, double threshold, const vector<unsigned int>& reference, size_t nbSequencesRef)
{
  AlignmentPartitionScores2& nodeScores = scores[node.getId()]; //This was created during the first pass.
  size_t nbSites = sites.getNumberOfSites();
  nodeScores.up.assign(nbSites, 0);
  nodeScores.nbSitesUp = 0;
  nodeScores.nbCharsUp = 0;
  if (node.hasFather()) {
    //If the node is root, there is nothing to do here...
    const Node* father = node.getFather();
    AlignmentPartitionScores2& fatherScores = scores[father->getId()];
    //We initialize the up array to the one of it father:
    nodeScores.up = fatherScores.up;
    nodeScores.nbSequencesUp = fatherScores.nbSequencesUp;
    //No we check all uncle nodes:
    for (size_t k = 0; k < father->getNumberOfSons(); ++k) {
      const Node* uncle = father->getSon(k);
      if (uncle != &node) {
        AlignmentPartitionScores2& uncleScores = scores[uncle->getId()];
        nodeScores.nbSequencesUp += uncleScores.nbSequencesDown;
        for (size_t i = 0; i < nbSites; ++i) {
          nodeScores.up[i] += uncleScores.down[i];
        }
      }
    }
  } else {
    nodeScores.nbSequencesUp = 0;
  }
  //Finally we count all sites:
  for (size_t i = 0; i < nbSites; ++i) {
    if (static_cast<double>(nodeScores.up[i] + reference[i]) / static_cast<double>(nodeScores.nbSequencesUp + nbSequencesRef) >= threshold) {
      nodeScores.nbSitesUp++;
      nodeScores.nbCharsUp += nodeScores.up[i];
    }
  }

  // Recursive call on son nodes (this does not do anything for leaves):
  for (size_t k = 0; k < node.getNumberOfSons(); ++k) {
    gapOptimizerSecondTreeTraversal2(*node.getSon(k), sites, scores, threshold, reference, nbSequencesRef);
  }
}

// Update function.
void gapOptimizerUpstreamUpdate2(const Node& node, const Node* from, const SiteContainer& sites, Index2& scores, double threshold, const vector<unsigned int>& reference, size_t nbSequencesRef)
{
  if (!node.isLeaf()) {
    size_t nbSites = sites.getNumberOfSites();
    AlignmentPartitionScores2& nodeScores = scores[node.getId()];
    nodeScores.down.assign(nbSites, 0);
    nodeScores.nbSequencesDown = 0;
    nodeScores.nbSitesDown = 0;
    nodeScores.nbCharsDown = 0;
    for (size_t k = 0; k < node.getNumberOfSons(); ++k) {
      //We loop over all son nodes:
      const Node* son = node.getSon(k);
      //We need to make computations for the subtree first and get the results:
      AlignmentPartitionScores2& sonScores = scores[son->getId()];
      //Now we make computations for this node:
      nodeScores.nbSequencesDown += sonScores.nbSequencesDown;
      for (size_t i = 0; i < nbSites; ++i)
        nodeScores.down[i] += sonScores.down[i];
      if (from && son->getId() != from->getId()) {
        //Now we update all up nodes in uncle nodes:
        gapOptimizerSecondTreeTraversal2(*son, sites, scores, threshold, reference, nbSequencesRef);
      }
    }
    //Finally we count all sites:
    for (size_t i = 0; i < nbSites; ++i) {
      if (static_cast<double>(nodeScores.down[i] + reference[i]) / static_cast<double>(nodeScores.nbSequencesDown + nbSequencesRef) >= threshold) {
        nodeScores.nbSitesDown++;
        nodeScores.nbCharsDown += nodeScores.down[i];
      }
    }
  }
  if (node.hasFather()) gapOptimizerUpstreamUpdate2(*node.getFather(), &node, sites, scores, threshold, reference, nbSequencesRef);
}

AlignmentPartitionScores2& gapOptimizerDownstreamUpdate2(const Node& node, const SiteContainer& sites, Index2& scores, double threshold, bool filterGaps, bool filterUnresolved, const vector<unsigned int>& reference, size_t nbSequencesRef)
{
  AlignmentPartitionScores2& nodeScores = scores[node.getId()]; //This eventually creates scores for this node.
  size_t nbSites = sites.getNumberOfSites();
  nodeScores.id = node.getId();
  nodeScores.down.assign(nbSites, 0);
  nodeScores.nbSitesDown = 0;
  nodeScores.nbCharsDown = 0;
  
  //We distinguish two cases, whether the node is a leaf or an inner node:
  if (node.isLeaf()) {
    //This is a leaf, we initialize arrays from the sequence data.
    const Sequence& seq = sites.getSequence(node.getName()); //We assume that the alignment contains all leaves of the tree.
    //This one at least is trivial:
    nodeScores.nbSequencesDown = 1;
    //Now we need to initialize the bit vector and count sites:
    for (size_t i = 0; i < nbSites; ++i) {
      bool test = !((filterGaps && seq.getAlphabet()->isGap(seq[i])) || (filterUnresolved && seq.getAlphabet()->isUnresolved(seq[i])));
      nodeScores.down[i] = (test ? 1 : 0);
    }
  } else {
    //This is an inner node.
    //All calculations are performed recursively from the results of son nodes:
    nodeScores.nbSequencesDown = 0;
    for (size_t k = 0; k < node.getNumberOfSons(); ++k) {
      //We loop over all son nodes:
      const Node* son = node.getSon(k);
      //We need to make computations for the subtree first and get the results:
      AlignmentPartitionScores2& sonScores = scores[son->getId()];
      //Now we make computations for this node:
      nodeScores.nbSequencesDown += sonScores.nbSequencesDown;
      for (size_t i = 0; i < nbSites; ++i)
        nodeScores.down[i] += sonScores.down[i];
    }
  }
  //Finally we count all sites:
  for (size_t i = 0; i < nbSites; ++i) {
    if (static_cast<double>(nodeScores.down[i] + reference[i]) / static_cast<double>(nodeScores.nbSequencesDown + nbSequencesRef) >= threshold) {
      nodeScores.nbSitesDown++;
      nodeScores.nbCharsDown += nodeScores.down[i];
    }
  }
  return nodeScores;
}



class Selector {
  public:
    virtual size_t getSelection(size_t nbChoices) const = 0;
};

class AutoSelector: public Selector {
  public:
    virtual size_t getSelection(size_t nbChoices) const { return (nbChoices > 0 ? 1 : 0); }
};

class InputSelector: public Selector {
  public:
    virtual size_t getSelection(size_t nbChoices) const {
      cout << "Your choice (0 to stop):" << endl;
      size_t choice;
      cin >> choice;
      return choice;
    }
};



int main(int args, char ** argv)
{
  cout << "******************************************************************" << endl;
  cout << "*           Bio++ Alignment Optimizer, version 1.0.0.            *" << endl;
  cout << "* Author: J. Dutheil                        Last Modif. 26/05/17 *" << endl;
  cout << "*         E. Figuet                                              *" << endl;
  cout << "******************************************************************" << endl;
  cout << endl;
  
  if (args == 1)
  {
    help();
    return 0;
  }
  
  try {

  BppApplication bppalnoptim(args, argv, "BppAlnOpt");
  bppalnoptim.startTimer();

  //Get sequences:
  auto_ptr<Alphabet> alphabet(SequenceApplicationTools::getAlphabet(bppalnoptim.getParams()));
  auto_ptr<SiteContainer> sites(SequenceApplicationTools::getSiteContainer(alphabet.get(), bppalnoptim.getParams()));

  //Get options:
  double threshold = ApplicationTools::getDoubleParameter("threshold", bppalnoptim.getParams(), 0.5, "", true, 1);
  ApplicationTools::displayResult("Minimum amount of data per site", threshold);
  bool filterGaps = ApplicationTools::getBooleanParameter("filter_gaps", bppalnoptim.getParams(), true, "", true, 1);
  ApplicationTools::displayBooleanResult("Filter gap characters", filterGaps);
  bool filterUnresolved = ApplicationTools::getBooleanParameter("filter_unresolved", bppalnoptim.getParams(), false, "", true, 1);
  ApplicationTools::displayBooleanResult("Filter unresolved characters", filterUnresolved);
  if (!filterGaps && !filterUnresolved)
    throw Exception("Error, either gap, unresolved characters or both should be filtered!");

  //Deal with reference alignment, if any.
  vector<string> refSequencesNames = ApplicationTools::getVectorParameter<string>("reference.sequences", bppalnoptim.getParams(), ',', "", "", true, 1);
  size_t nbSequencesRef = refSequencesNames.size();
  ApplicationTools::displayResult("Number of reference sequences", nbSequencesRef);
  AlignedSequenceContainer refAln(alphabet.get());
  for (size_t i = 0; i < nbSequencesRef; ++i) {
    try {
      auto_ptr<Sequence> seq(sites->removeSequence(refSequencesNames[i]));
      refAln.addSequence(*seq);
    } catch(SequenceNotFoundException) {
      throw Exception("No sequence with name '" + refSequencesNames[i] + "' was found in the input alignment.");
    }
  }
  size_t totalNbSites = sites->getNumberOfSites();
  ApplicationTools::displayResult("Total number of sequences", sites->getNumberOfSequences());
  ApplicationTools::displayResult("Total number of sites", totalNbSites);
  vector<unsigned int> reference(totalNbSites, 0);
  for (size_t i = 0; i < nbSequencesRef; ++i) {
    const Sequence& seq = refAln.getSequence(i);
    for (size_t j = 0; j < totalNbSites; ++j) {
      if (!((filterGaps && alphabet->isGap(seq[j])) || (filterUnresolved && alphabet->isUnresolved(seq[j]))))
        reference[j]++;
    }
  }

  //Get tree:
  auto_ptr< TreeTemplate<Node> > tree;
  auto_ptr< TreeTemplate<Node> > origTree;
  
  //input or cluster:
  string inputTree = ApplicationTools::getStringParameter("input.tree.method", bppalnoptim.getParams(), "AutoCluster", "", false, 1);
  map<string, string> inputTreeParams;
  string inputTreeName;
  KeyvalTools::parseProcedure(inputTree, inputTreeName, inputTreeParams);
  if (inputTreeName == "AutoCluster") {
    //Compute pairwise distance matrix:
    size_t n = sites->getNumberOfSequences();
    size_t ns = sites->getNumberOfSites();

    //First we compress sequences:
    vector< vector<bool> > cseqs(n);
    ApplicationTools::displayTask("Compressing sequences", true);
    for (size_t i = 0; i < n; ++i) {
      ApplicationTools::displayGauge(i, n - 1);
      const Sequence& seq = sites->getSequence(i);
      cseqs[i].resize(ns);
      for (size_t j = 0; j < ns; ++j)
        cseqs[i][j] = !alphabet->isGap(seq[j]);
    }
    ApplicationTools::displayTaskDone();

    size_t totmem = (sizeof(vector< vector<double> >) + (sizeof(vector<double>) + sizeof(double) * n) * n) / 1024 / 1024;
    ApplicationTools::displayResult("Memory required to store distances (Mb)", totmem);
    ApplicationTools::displayTask("Computing pairwise overlap matrix", true);
    DistanceMatrix d(sites->getSequencesNames());
    size_t k = 0, m = n * (n - 1) / 2;
    for (size_t i = 1; i < n; ++i) {
      d(i, i) = 0.;
      for (size_t j = 0; j < i; ++j) {
        k++;
        ApplicationTools::displayGauge(k, m);
        d(i, j) = d(j, i) = overlapDistance(cseqs[i], cseqs[j]);
      }
    }
    ApplicationTools::displayTaskDone();
    //Free mem:
    cseqs.clear();

    //Now we build a cluster tree:
    string linkage = ApplicationTools::getStringParameter("linkage", inputTreeParams, "median", "", true, 2);
    string linkMode = "";
    if (linkage == "complete") linkMode = HierarchicalClustering::COMPLETE; 
    if (linkage == "single"  ) linkMode = HierarchicalClustering::SINGLE; 
    if (linkage == "average" ) linkMode = HierarchicalClustering::AVERAGE;
    if (linkage == "median"  ) linkMode = HierarchicalClustering::MEDIAN;
    if (linkage == "ward"    ) linkMode = HierarchicalClustering::WARD;
    if (linkage == "centroid") linkMode = HierarchicalClustering::CENTROID;
    if (linkMode == "") throw Exception("Invalid linkage mode: " + linkage + ".");
    ApplicationTools::displayResult("Clustering linkage mode", linkage);
    ApplicationTools::displayTask("Computing cluster tree", true);
    HierarchicalClustering hclust(linkMode, d, true);
    tree.reset(hclust.getTree());
    tree->unroot();
    ApplicationTools::displayTaskDone();
  } else if (inputTreeName == "File") {
    tree.reset(dynamic_cast<TreeTemplate<Node> *>(PhylogeneticsApplicationTools::getTree(bppalnoptim.getParams())));
    if (tree->isRooted()) {
      ApplicationTools::displayWarning("Tree has been unrooted.");
      tree->unroot();
    }
    origTree.reset(tree->clone()); //We keep a copy of the original tree
    for (size_t i = 0; i < nbSequencesRef; ++i) {
      try  {
        TreeTemplateTools::dropLeaf(*tree, refSequencesNames[i]);
      } catch(NodeNotFoundException&) {} //Just ignore this, tree does not include ref sequences and that is ok.
    }
  } else {
    throw Exception("Unrecognized tree input method: " + inputTreeName + ".");
  }
  //Reorder the alignment according to tree and use sequence storage instead of site storage (more efficient);
  vector<string> leafNames = tree->getLeavesNames();
  AlignedSequenceContainer* asc = new AlignedSequenceContainer(alphabet.get());
  for (size_t i = 0; i < leafNames.size(); ++i) {
    asc->addSequence(sites->getSequence(leafNames[i]));
  }
  sites.reset(asc);

  //The main loop:
  string logPath = ApplicationTools::getAFilePath("output.log", bppalnoptim.getParams(), false, false, "", true, "bppalnoptim.log", 1);
  ofstream logfile(logPath.c_str());  
  logfile << "Iteration\tChoice\tNode\tNbSequences\tDiffSequences\tNbSites\tDiffSites\tNbChars\tDiffChars\tMeanEntropy" << endl;

  //Build initial scores. This will be updated after each iteration, if needed.
  Index2 indexedScores;
  gapOptimizerFirstTreeTraversal2(*tree->getRootNode(), *sites, indexedScores, threshold, filterGaps, filterUnresolved, reference, nbSequencesRef); 
  gapOptimizerSecondTreeTraversal2(*tree->getRootNode(), *sites, indexedScores, threshold, reference, nbSequencesRef);

  //Choose analysis mode:
  string methodDesc = ApplicationTools::getStringParameter("method", bppalnoptim.getParams(), "Input(show=10)", "", false, 1);
  string methodName;
  map<string, string> methodArgs;
  KeyvalTools::parseProcedure(methodDesc, methodName, methodArgs);

  //We try to maximize the number of sites while removing the minimum number of sequences
  size_t nbDisplay = 1;
  auto_ptr<Selector> selector;
  unsigned int minNbSequences = 0;
  unsigned int minNbSitesRequired = static_cast<unsigned int>(sites->getNumberOfSites()) + 1;
  if (methodName == "Auto" || methodName == "Diagnostic") {
    selector.reset(new AutoSelector());
    if (methodName == "Auto") {
      minNbSequences = ApplicationTools::getParameter<unsigned int>("min_nb_sequences", methodArgs, 0, "", true, 1);
      if (minNbSequences == 0) {
        double f = ApplicationTools::getParameter<double>("min_relative_nb_sequences", methodArgs, 0, "", true, 1);
        minNbSequences = static_cast<unsigned int>(floor(static_cast<double>(nbSequencesRef + sites->getNumberOfSequences()) * f));
      }
      minNbSitesRequired = ApplicationTools::getParameter<unsigned int>("min_nb_sites", methodArgs, minNbSitesRequired, "", true, 1);
      if (minNbSitesRequired == 0) {
        double f = ApplicationTools::getParameter<double>("min_relative_nb_sites", methodArgs, 1, "", true, 1);
        minNbSitesRequired = static_cast<unsigned int>(floor(static_cast<double>(sites->getNumberOfSites()) * f));
      }
    }
  } else if (methodName == "Input") {
    selector.reset(new InputSelector());
    nbDisplay = ApplicationTools::getParameter<size_t>("show", methodArgs, 10, "", false, 2);
  } else {
    throw Exception("Unrecognized selection method: " + methodDesc);
  }
  if (minNbSequences > 0)
    ApplicationTools::displayResult("Minimum number of sequences to keep", minNbSequences);
  if (minNbSitesRequired <= static_cast<unsigned int>(sites->getNumberOfSites()))
    ApplicationTools::displayResult("Stop when this number of sites is reached", minNbSitesRequired);

  auto_ptr<Comparator> comp;
  string compCrit = ApplicationTools::getStringParameter("comparison", bppalnoptim.getParams(), "MaxSites", "", false, 1);
  if (compCrit == "MaxSites") {
    comp.reset(new Comparator('a'));
  } else if (compCrit == "MaxChars") {
    comp.reset(new Comparator('b'));
  } else {
    throw Exception("Unrecognized comparison criterion: " + compCrit);
  }
  ApplicationTools::displayResult("Comparison criterion:", compCrit);

  bool test = true;
  unsigned int iteration = 0;
  vector<string> removedSequenceNames;
  double meanE = 0;
  for (size_t i = 0; i < sites->getNumberOfSites(); ++i) {
    meanE += SiteTools::variabilityShannon(sites->getSite(i), false);
  }
  meanE /= static_cast<double>(sites->getNumberOfSites());
  logfile << "0\t0\t" << tree->getRootId() << "\t" << indexedScores[tree->getRootId()].nbSequencesDown << "\t0\t" << indexedScores[tree->getRootId()].nbSitesDown << "\t0\t" << indexedScores[tree->getRootId()].nbCharsDown << "\t0\t" << meanE << endl;
  while (test) {
    iteration++;
    unsigned int nbSequences = indexedScores[tree->getRootId()].nbSequencesDown;
    unsigned int nbSites = indexedScores[tree->getRootId()].nbSitesDown;
    unsigned int nbChars = indexedScores[tree->getRootId()].nbCharsDown;
    ApplicationTools::displayResult("Number of sequences in alignment", nbSequences);
    ApplicationTools::displayResult("Number of sites in alignment", nbSites);
    ApplicationTools::displayResult("Number of chars in alignment", nbChars);
   
    //Check if we have enough sites:
    if (nbSites >= minNbSitesRequired) {
      test = false;
      continue;
    }
    ApplicationTools::displayResult("Mean site entropy in alignment", meanE);

    //Maximizes the number of sites:
    //First get all splits which improve the number of sites and meet the stopping condition, if any:
    vector<Group> improvingGroups;
    for (Index2::iterator it = indexedScores.begin(); it != indexedScores.end(); ++it) {
      Group gup = it->second.getGroupUp();
      if (comp->isImproving(gup, nbSites, nbChars) && gup.nbSequences >= minNbSequences)
        improvingGroups.push_back(gup);
      Group gdo = it->second.getGroupDown(); 
      if (comp->isImproving(gdo, nbSites, nbChars) && gdo.nbSequences >= minNbSequences)
        improvingGroups.push_back(gdo); 
    }
    //Second sort according to number of sequences removed:
    sort(improvingGroups.begin(), improvingGroups.end(), *comp.get());

    //Show all beneficial splits:

    cout << "Choice\tNode\t#seq\t%seq\t#sites\t%sites\t#chars\t%chars" << endl;
    size_t nbChoices = nbSequences > 2 ? min(improvingGroups.size(), nbDisplay) : 0;
    for (size_t i = 0; i < nbChoices; ++i) {
      Group& group = improvingGroups[i];
      cout << (i+1) << ")\t" << group.id << "\t" << group.nbSequences << "\t" << (static_cast<double>(group.nbSequences) - static_cast<double>(nbSequences)) << "\t" << group.nbSites << "\t+" << (group.nbSites - nbSites) << "\t" << group.nbChars << "\t" << (static_cast<double>(group.nbChars) - static_cast<double>(nbChars)) << endl;
    }
    size_t choice = selector->getSelection(nbChoices);

    if (choice == 0)
      test = false;
    else {
      Group& group = improvingGroups[choice - 1];
      Node* node = tree->getNode(group.id);
      Node* father = node->getFather();
      if (group.direction == 'u') {
        //Main group is 'Up', 'Down' should be removed...
        Node* grandFather = 0;
        if (father->hasFather()) 
          grandFather = father->getFather();
        vector<string> removedLeaves = TreeTemplateTools::getLeavesNames(*node);
        TreeTemplateTools::dropSubtree(*tree, node);
        removedSequenceNames.insert(removedSequenceNames.end(), removedLeaves.begin(), removedLeaves.end());
        VectorSiteContainer* selectedSites = new VectorSiteContainer(sites->getAlphabet());
        SequenceContainerTools::getSelectedSequences(*sites, tree->getLeavesNames(), *selectedSites);
        sites.reset(selectedSites);
        //We need to recompute down scores of all parent nodes up to the root:
        if (grandFather) {
          gapOptimizerUpstreamUpdate2(*grandFather, 0, *sites, indexedScores, threshold, reference, nbSequencesRef);
          gapOptimizerSecondTreeTraversal2(*grandFather, *sites, indexedScores, threshold, reference, nbSequencesRef);
        } else {
          //The remove clade was branching at the root
          if (tree->isRooted()) {
            tree->unroot(); //The tree may have become artifically rooted as we removed one son node, there may be only two left :(
            gapOptimizerDownstreamUpdate2(*tree->getRootNode(), *sites, indexedScores, threshold, filterGaps, filterUnresolved, reference, nbSequencesRef);
          }
          gapOptimizerSecondTreeTraversal2(*tree->getRootNode(), *sites, indexedScores, threshold, reference, nbSequencesRef);
        }
      } else {
        //Main group is 'Down', 'Up' should be removed...
        //We set the node as new root node and delete the rest of the tree.
        father->removeSon(node);
        vector<string> removedLeaves = TreeTemplateTools::getLeavesNames(*tree->getRootNode());
        TreeTemplateTools::deleteSubtree(tree->getRootNode());
        removedSequenceNames.insert(removedSequenceNames.end(), removedLeaves.begin(), removedLeaves.end());
        tree->setRootNode(node);
        if (tree->getNumberOfLeaves() > 1) {
          tree->unroot();
          node = tree->getRootNode(); //update with actual root node
          //We need to recompute 'down' scores at this node:
          AlignmentPartitionScores2& nodeScores = indexedScores[tree->getRootId()];
          size_t alnSize = sites->getNumberOfSites();
          nodeScores.down.assign(alnSize, 0);
          nodeScores.nbSitesDown = 0;
          nodeScores.nbCharsDown = 0;
          nodeScores.nbSequencesDown = 0;
          for (size_t k = 0; k < node->getNumberOfSons(); ++k) {
            //We loop over all son nodes:
            const Node* son = node->getSon(k);
            AlignmentPartitionScores2& sonScores = indexedScores[son->getId()];
            //Now we make computations for this node:
            nodeScores.nbSequencesDown += sonScores.nbSequencesDown;
            for (size_t i = 0; i < alnSize; ++i)
              nodeScores.down[i] += sonScores.down[i];
            }
            //Finally we count all sites:
            for (size_t i = 0; i < alnSize; ++i) {
            if (static_cast<double>(nodeScores.down[i]) / static_cast<double>(nodeScores.nbSequencesDown) >= threshold) {
              nodeScores.nbSitesDown++;
              nodeScores.nbCharsDown += nodeScores.down[i];
            }
          }
        }
        VectorSiteContainer* selectedSites = new VectorSiteContainer(sites->getAlphabet());
        SequenceContainerTools::getSelectedSequences(*sites, tree->getLeavesNames(), *selectedSites);
        sites.reset(selectedSites);
        //We recompute up scores only:
        gapOptimizerSecondTreeTraversal2(*tree->getRootNode(), *sites, indexedScores, threshold, reference, nbSequencesRef);
      }
      //We update the index by removing unecessary nodes:
      vector<int> ids = tree->getNodesId();
      for (Index2::iterator it = indexedScores.begin(); it != indexedScores.end();) {
        if (find(ids.begin(), ids.end(), it->first) == ids.end()) {
          indexedScores.erase(it++);
        } else {
          ++it;
        }
      }
   
      //Compute entropy:
      meanE = 0;
      for (size_t i = 0; i < sites->getNumberOfSites(); ++i) {
        meanE += SiteTools::variabilityShannon(sites->getSite(i), false);
      }
      meanE /= static_cast<double>(sites->getNumberOfSites()); 
      
      //Write to log file:
      logfile << iteration << "\t" << choice << "\t" << group.id << "\t" << group.nbSequences << "\t" << (static_cast<double>(group.nbSequences) - static_cast<double>(nbSequences)) << "\t" << group.nbSites << "\t+" << (group.nbSites - nbSites) << "\t" << group.nbChars << "\t" << (static_cast<double>(group.nbChars) - static_cast<double>(nbChars)) << "\t" << meanE << endl;
    }
  }
  logfile.close();

  if (methodName != "Diagnostic") {
    //Write final alignment:
    if (ApplicationTools::getStringParameter("output.sequence.file", bppalnoptim.getParams(), "none", "", true, 1) != "none") {
      if (nbSequencesRef > 0) {
        SequenceContainerTools::append(refAln, *sites.get());
        SequenceApplicationTools::writeAlignmentFile(refAln, bppalnoptim.getParams());
      } else {
        SequenceApplicationTools::writeAlignmentFile(*sites.get(), bppalnoptim.getParams());
      }
    }

    //Write final tree:
    if (origTree.get()) {
      for (size_t i = 0; i < removedSequenceNames.size(); ++i)
        TreeTemplateTools::dropLeaf(*origTree, removedSequenceNames[i]);
      PhylogeneticsApplicationTools::writeTree(*origTree.get(), bppalnoptim.getParams());
    } else {
      //Warning! This will not contain the reference sequence names! (and that is ok :) ).
      PhylogeneticsApplicationTools::writeTree(*tree.get(), bppalnoptim.getParams());
    }
  }

  bppalnoptim.done();
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

