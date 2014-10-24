//
// File: bppAlignmentOptimizer.cpp
// Created by: Julien Dutheil
// Created on: Sunday, December 2nd 2007 16:48
//

/*
Copyright or © or Copr. Bio++ Development Team

This software is a computer program whose purpose is to estimate
phylogenies and evolutionary parameters from a dataset according to
the maximum likelihood principle.

This software is governed by the CeCILL  license under French law and
abiding by the rules of distribution of free software.  You can  use, 
modify and/ or redistribute the software under the terms of the CeCILL
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info". 

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability. 

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or 
data to be ensured and,  more generally, to use and operate it in the 
same conditions as regards security. 

The fact that you are presently reading this means that you have had
knowledge of the CeCILL license and that you accept its terms.
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
#include <Bpp/Seq/Alphabet.all>
#include <Bpp/Seq/Container.all>
#include <Bpp/Seq/App/SequenceApplicationTools.h>

// From bpp-phyl:
#include <Bpp/Phyl/TreeTemplate.h>
#include <Bpp/Phyl/App/PhylogeneticsApplicationTools.h>

using namespace bpp;

void help()
{
  (*ApplicationTools::message << "__________________________________________________________________________").endLine();
  (*ApplicationTools::message << "bppalnoptim parameter1_name=parameter1_value").endLine();
  (*ApplicationTools::message << "      parameter2_name=parameter2_value ... param=option_file").endLine();
  (*ApplicationTools::message).endLine();
  (*ApplicationTools::message << "  Refer to the Bio++ Program Suite Manual for a list of available options.").endLine();
  (*ApplicationTools::message << "__________________________________________________________________________").endLine();
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
      nbSitesUp(0),
      nbSitesDown(0),
      nbSequencesUp(0),
      nbSequencesDown(0),
      nbCharsUp(0),
      nbCharsDown(0)
    {}

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


// Recursive function.
// This will initialize the scores map and set the down variables.
// For convenience, we return the scores for the root subtree.
AlignmentPartitionScores1& gapOptimizerFirstTreeTraversal1(const Node& node, const SiteContainer& sites, map<int, AlignmentPartitionScores1>& scores)
{
  AlignmentPartitionScores1& nodeScores = scores[node.getId()]; //This creates scores for this node.
  size_t nbSites = sites.getNumberOfSites();
  nodeScores.id = node.getId();
  nodeScores.up.assign(nbSites, true);
  nodeScores.down.assign(nbSites, true);
  
  //We distinguish two cases, whether the node is a leaf or an inner node:
  if (node.isLeaf()) {
    //This is a leaf, we initialize arrays from the sequence data.
    const Sequence& seq = sites.getSequence(node.getName()); //We assume that the alignment contains all leaves of the tree.
    //This one at least is trivial:
    nodeScores.nbSequencesDown = 1;
    //Now we need to initialize the bit vector and count sites:
    for (size_t i = 0; i < nbSites; ++i) {
      bool test = !seq.getAlphabet()->isGap(seq[i]);
      nodeScores.down[i] = test; //NB: we need to adapt in case we also want to ignore generic characters.
      if (test) nodeScores.nbSitesDown++;
    }
  } else {
    //This is an inner node.
    //All calculations are performed recursively from the results of son nodes:
    for (size_t k = 0; k < node.getNumberOfSons(); ++k) {
      //We loop over all son nodes:
      const Node* son = node.getSon(k);
      //We need to make computations for the subtree first and get the results:
      AlignmentPartitionScores1& sonScores = gapOptimizerFirstTreeTraversal1(*son, sites, scores);
      //Now we make computations for this node:
      nodeScores.nbSequencesDown += sonScores.nbSequencesDown;
      for (size_t i = 0; i < nbSites; ++i)
        nodeScores.down[i] = nodeScores.down[i] && sonScores.down[i];
    }
    //Finally we count all sites:
    for (size_t i = 0; i < nbSites; ++i)
      if (nodeScores.down[i]) nodeScores.nbSitesDown++;
    nodeScores.nbCharsDown = nodeScores.nbSequencesDown * nodeScores.nbSitesDown;
  }
  return nodeScores;
}

// Recursive function.
// This will initialize the scores map and set the up variables.
// This must be ran after the first tree traversal, so that all scores are already initialized.
void gapOptimizerSecondTreeTraversal1(const Node& node, const SiteContainer& sites, map<int, AlignmentPartitionScores1>& scores)
{
  AlignmentPartitionScores1& nodeScores = scores[node.getId()]; //This was created during the first pass.
  size_t nbSites = sites.getNumberOfSites();
  if (node.hasFather()) {
    const Node* father = node.getFather();
    //If the node is root, there is nothing to do here...
    AlignmentPartitionScores1& fatherScores = scores[father->getId()];
    //We initialize the up array to the one of it father:
    nodeScores.up = fatherScores.up;
    nodeScores.nbSequencesUp = fatherScores.nbSequencesUp;
    //No we check all uncle nodes:
    for (size_t k = 0; k < father->getNumberOfSons(); ++k) {
      const Node* uncle = father->getSon(k);
      if (uncle != &node) {
        AlignmentPartitionScores1& uncleScores = scores[uncle->getId()];
        nodeScores.nbSequencesUp += uncleScores.nbSequencesDown;
        for (size_t i = 0; i < nbSites; ++i) {
          nodeScores.up[i] = nodeScores.up[i] && uncleScores.down[i];
        }
      }
    }
  } else {
    //This is the root node, the up vector equals the down one:
    nodeScores.up = nodeScores.down;
  }
  //Finally we count all sites:
  for (size_t i = 0; i < nbSites; ++i)
    if (nodeScores.up[i]) nodeScores.nbSitesUp++;
  nodeScores.nbCharsUp = nodeScores.nbSequencesUp * nodeScores.nbSitesUp;

  // Recursive call on son nodes (this does not do anything for leaves):
  for (size_t k = 0; k < node.getNumberOfSons(); ++k) {
    gapOptimizerSecondTreeTraversal1(*node.getSon(k), sites, scores);
  }
}

#define Index2 map<int, AlignmentPartitionScores2>

// Recursive function.
// This will initialize the scores map and set the down variables.
// For convenience, we return the scores for the root subtree.
AlignmentPartitionScores2& gapOptimizerFirstTreeTraversal2(const Node& node, const SiteContainer& sites, Index2& scores, double threshold, bool filterUnresolved)
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
      bool test = !seq.getAlphabet()->isGap(seq[i]);
      if (filterUnresolved)
        test = test & !seq.getAlphabet()->isUnresolved(seq[i]);
      nodeScores.down[i] = (test ? 1 : 0); //NB: we need to adapt in case we also want to ignore generic characters.
      if (test) nodeScores.nbSitesDown++;
    }
  } else {
    //This is an inner node.
    //All calculations are performed recursively from the results of son nodes:
    nodeScores.nbSequencesDown = 0;
    for (size_t k = 0; k < node.getNumberOfSons(); ++k) {
      //We loop over all son nodes:
      const Node* son = node.getSon(k);
      //We need to make computations for the subtree first and get the results:
      AlignmentPartitionScores2& sonScores = gapOptimizerFirstTreeTraversal2(*son, sites, scores, threshold, filterUnresolved);
      //Now we make computations for this node:
      nodeScores.nbSequencesDown += sonScores.nbSequencesDown;
      for (size_t i = 0; i < nbSites; ++i)
        nodeScores.down[i] += sonScores.down[i];
    }
    //Finally we count all sites:
    for (size_t i = 0; i < nbSites; ++i) {
      if (static_cast<double>(nodeScores.down[i]) / static_cast<double>(nodeScores.nbSequencesDown) >= threshold) {
        nodeScores.nbSitesDown++;
        nodeScores.nbCharsDown += nodeScores.down[i];
      }
    }
  }
  return nodeScores;
}

// Recursive function.
// This will initialize the scores map and set the up variables.
// This must be ran after the first tree traversal, so that all scores are already initialized.
void gapOptimizerSecondTreeTraversal2(const Node& node, const SiteContainer& sites, Index2& scores, double threshold)
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
    if (static_cast<double>(nodeScores.up[i]) / static_cast<double>(nodeScores.nbSequencesUp) >= threshold) {
      nodeScores.nbSitesUp++;
      nodeScores.nbCharsUp += nodeScores.up[i];
    }
  }

  // Recursive call on son nodes (this does not do anything for leaves):
  for (size_t k = 0; k < node.getNumberOfSons(); ++k) {
    gapOptimizerSecondTreeTraversal2(*node.getSon(k), sites, scores, threshold);
  }
}

// Update function.
void gapOptimizerUpstreamUpdate2(const Node& node, const Node* from, const SiteContainer& sites, Index2& scores, double threshold)
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
        gapOptimizerSecondTreeTraversal2(*son, sites, scores, threshold);
      }
    }
    //Finally we count all sites:
    for (size_t i = 0; i < nbSites; ++i) {
      if (static_cast<double>(nodeScores.down[i]) / static_cast<double>(nodeScores.nbSequencesDown) >= threshold) {
        nodeScores.nbSitesDown++;
        nodeScores.nbCharsDown += nodeScores.down[i];
      }
    }
  }
  if (node.hasFather()) gapOptimizerUpstreamUpdate2(*node.getFather(), &node, sites, scores, threshold);
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
  cout << "*           Bio++ Alignment Optimizer, version 2.2.0.            *" << endl;
  cout << "* Author: J. Dutheil                        Last Modif. 02/10/14 *" << endl;
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

  //Get tree:
  auto_ptr< TreeTemplate<Node> > tree;
  tree.reset(dynamic_cast<TreeTemplate<Node> *>(PhylogeneticsApplicationTools::getTree(bppalnoptim.getParams())));

  //Get options:
  double threshold = ApplicationTools::getDoubleParameter("threshold", bppalnoptim.getParams(), 0.5, "", true, 1);
  ApplicationTools::displayResult("Minimum amount of data per site", threshold);
  bool filterUnresolved = ApplicationTools::getBooleanParameter("filter_unresolved", bppalnoptim.getParams(), false, "", true, 1);
  ApplicationTools::displayBooleanResult("Filter unresolved characters", filterUnresolved);

  //The main loop:
  string logPath = ApplicationTools::getAFilePath("output.log", bppalnoptim.getParams(), false, false, "", true, "bppalnoptim.log", 1);
  ofstream logfile(logPath.c_str());  
  logfile << "Iteration\tChoice\tNode\tNbSequences\tDiffSequences\tNbSites\tDiffSites\tNbChars\tDiffChars" << endl;

  //Build initial scores. This will be updated after each iteration, if needed.
  Index2 indexedScores;
  gapOptimizerFirstTreeTraversal2(*tree->getRootNode(), *sites, indexedScores, threshold, filterUnresolved); 
  gapOptimizerSecondTreeTraversal2(*tree->getRootNode(), *sites, indexedScores, threshold);

  //Choose analysis mode:
  string methodDesc = ApplicationTools::getStringParameter("method", bppalnoptim.getParams(), "Input(show=10)", "", false, 1);
  string methodName;
  map<string, string> methodArgs;
  KeyvalTools::parseProcedure(methodDesc, methodName, methodArgs);

  //We try to maximize the number of sites while removing the minimum number of sequences
  size_t nbDisplay = 1;
  auto_ptr<Selector> selector;
  if (methodName == "Auto" || methodName == "Diagnostic") {
    selector.reset(new AutoSelector());
  } else if (methodName == "Input") {
    selector.reset(new InputSelector());
    nbDisplay = ApplicationTools::getParameter<size_t>("show", methodArgs, 10, "", false, 2);
  } else {
    throw Exception("Unrecognized selection method: " + methodDesc);
  }

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
  while (test) {
    iteration++;
    unsigned int nbSequences = indexedScores[tree->getRootId()].nbSequencesDown;
    unsigned int nbSites = indexedScores[tree->getRootId()].nbSitesDown;
    unsigned int nbChars = indexedScores[tree->getRootId()].nbCharsDown;
    ApplicationTools::displayResult("Number of sequences in alignment", nbSequences);
    ApplicationTools::displayResult("Number of sites in alignment", nbSites);
    ApplicationTools::displayResult("Number of chars in alignment", nbChars);

    //Maximizes the number of sites:
    //First ge all splits which improve the number of sites:
    vector<Group> improvingGroups;
    for (Index2::iterator it = indexedScores.begin(); it != indexedScores.end(); ++it) {
      Group gup = it->second.getGroupUp();
      if (comp->isImproving(gup, nbSites, nbChars))
        improvingGroups.push_back(gup);
      Group gdo = it->second.getGroupDown(); 
      if (comp->isImproving(gdo, nbSites, nbChars))
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
      logfile << iteration << "\t" << choice << "\t" << group.id << "\t" << group.nbSequences << "\t" << (static_cast<double>(group.nbSequences) - static_cast<double>(nbSequences)) << "\t" << group.nbSites << "\t+" << (group.nbSites - nbSites) << "\t" << group.nbChars << "\t" << (static_cast<double>(group.nbChars) - static_cast<double>(nbChars)) << endl;
      if (group.direction == 'u') {
        //Main group is 'Up', 'Down' should be removed...
        Node* grandFather = 0;
        if (father->hasFather()) 
          grandFather = father->getFather();
        TreeTemplateTools::dropSubtree(*tree, node);
        VectorSiteContainer* selectedSites = new VectorSiteContainer(sites->getAlphabet());
        SequenceContainerTools::getSelectedSequences(*sites, tree->getLeavesNames(), *selectedSites);
        sites.reset(selectedSites);
        //We need to recompute down scores of all parent nodes up to the root:
        if (grandFather) {
          gapOptimizerUpstreamUpdate2(*grandFather, 0, *sites, indexedScores, threshold);
          gapOptimizerSecondTreeTraversal2(*grandFather, *sites, indexedScores, threshold);
        } else {
          //The remove clade was branching at the root
          gapOptimizerSecondTreeTraversal2(*tree->getRootNode(), *sites, indexedScores, threshold);
        }
      } else {
        //Main group is 'Down', 'Up' should be removed...
        //We set the node as new root node and delete the rest of the tree.
        father->removeSon(node);
        TreeTemplateTools::deleteSubtree(tree->getRootNode());
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
        gapOptimizerSecondTreeTraversal2(*tree->getRootNode(), *sites, indexedScores, threshold);
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
    }
  }
  logfile.close();

  if (methodName != "Diagnostic") {
    //Write final alignment:
    SequenceApplicationTools::writeAlignmentFile(*sites.get(), bppalnoptim.getParams());

    //Write final tree:
    PhylogeneticsApplicationTools::writeTree(*tree.get(), bppalnoptim.getParams());
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

