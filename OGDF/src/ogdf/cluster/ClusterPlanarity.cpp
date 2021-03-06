/** \file
 * \brief Implementation of class ClusterPlanarity.
 *
 * \author Karsten Klein, Markus Chimani
 *
 * \par License:
 * This file is part of the Open Graph Drawing Framework (OGDF).
 * Copyright (C) 2005-2010
 *
 * \par
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * Version 2 or 3 as published by the Free Software Foundation;
 * see the file LICENSE.txt included in the packaging of this file
 * for details.
 *
 * \par
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * \par
 * You should have received a copy of the GNU General Public
 * License along with this program; if not, write to the Free
 * Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 * Boston, MA 02110-1301, USA.
 *
 * \see  http://www.gnu.org/copyleft/gpl.html
 ***************************************************************/

#include <ogdf/basic/basic.h>

#ifdef USE_ABACUS

#include <ogdf/cluster/ClusterPlanarity.h>
#include <ogdf/cluster/CconnectClusterPlanar.h>
#include <ogdf/basic/simple_graph_alg.h>
#include <ogdf/basic/extended_graph_alg.h>
#include <sstream>
//#include <ogdf/basic/exceptions.h>
#include <ogdf/internal/cluster/CPlanarity_Master.h>
#include <ogdf/planarity/BoyerMyrvold.h>

//#define writefeasiblegraphs
#ifdef OGDF_DEBUG
#include <ogdf/cluster/ClusterGraphAttributes.h>
#include <ogdf/fileformats/GraphIO.h>
#endif

using namespace abacus;

namespace ogdf {

struct connStruct {
	bool connected;
	node v1, v2;
	edge e;
};

bool ClusterPlanarity::isClusterPlanar(const ClusterGraph &CG) {
		List<nodePair> addedEdges;
		return isClusterPlanar(CG, addedEdges);
	}
bool ClusterPlanarity::isClusterPlanar(const ClusterGraph &CG, List<nodePair> &addedEdges) {
	m_optStatus = Master::Optimal;
	// We first check if there is more to do then just checking planarity on the
	// input graph.
	// Simple shortcut: With < 5 vertices, no non-planarity is possible...
	bool result = isPlanar(CG.constGraph());
	if (!result || (CG.numberOfClusters() == 1))
	{
		// Either non-planar or only root cluster exists, which does not restrict c-planarity.
		return result;
	}
	// We first create a copy of input G, and work solely on the copy

	// In case of the sm_new solution method, we partition the graph in
	// independent parts and test them separately
	// For all parts we test until non-c-planar or all tested.
	if (m_solmeth==sm_new)
	{
		// We use the ClusterAnalysis to search for independent bags
		// Here is the idea: We detect all bags that are minimum wrt
		// cluster inclusion (i.e. if a cluster contains a cluster c with
		// a single bag, we don't add the cluster c itself) but do
		// not contain an outeractive vertex wrt to smallest containing cluster
		// (i.e. the cluster used in the definition of bag).
		// The clustered subgraphs induced by these bags can be tested independently,
		// as we can move them freely in the drawing area of their enclosing parent cluster.

		ClusterAnalysis ca(CG, true); //Compute all structures, indyBags too
		// We can solve the c-planarity testing for all indyBags independently,
		// and in case all are c-planar, also our input c-graph is c-planar.
		const int numIndyBags = ca.numberOfIndyBags();
		GraphCopy** theGraphs = new GraphCopy * [numIndyBags]; //Stores copies for the bag graphs.
#ifdef OGDF_DEBUG
		cout << "Number of IndyBags "<<numIndyBags<<"\n";
#endif
		Logger::slout() << "Number of IndyBags "<<numIndyBags<<"\n";
		Array<List<node> > nodesInBag(numIndyBags); //Stores vertices for the bag graphs.
		const Graph & G = CG.constGraph();
		for(node v : G.nodes){
		   nodesInBag[ca.indyBagIndex(v)].pushBack(v);
		}

		for (int i = 0; i < numIndyBags && m_optStatus == Master::Optimal; i++)
		{
			// Create underlying graph
			theGraphs[i] = new GraphCopy();
			theGraphs[i]->createEmpty(G);
			// Judging from the interface and the description, there are two
			// methods in GraphCopy that allow to construct parts based on a
			// set of vertices, initByNodes and initByActiveNodes, where the
			// latter one seems to be appropriate and can be used with an
			// additional 3n work to initialize the NodeArray and mark the vertices.
			// However, even though the former is meant to be used for connected
			// components, it also works for set of connected components, and
			// an independent bag is such a creature.
			EdgeArray<edge> eCopy(G);
			theGraphs[i]->initByNodes(nodesInBag[i], eCopy);
			ClusterGraph bagCG(*theGraphs[i]);
			//node v;
			ClusterArray<List<node> > cNodes(CG);
			ClusterArray<List<cluster> > cChildren(CG);
			ClusterArray<cluster> cCopy(CG);
			// Run through all original vertices and store
			// lists of copies at each cluster that is part of the bag.
			// Note: We should not add an enclosing parent cluster below
			// root, i.e., when the root does only have a single child
			// and no vertices, we delete the child again.
			for(node u : nodesInBag[i])
			{
				cluster ct = CG.clusterOf(u);
				cNodes[ct].pushBack(theGraphs[i]->copy(u));
				// Check if we need to store the parent relation on the path
				// to the root. Indicator is: We have just added the first element.

				while ((ct != CG.rootCluster()) &&
						(cNodes[ct].size() + cChildren[ct].size() == 1))
				{
					cChildren[ct->parent()].pushBack(ct);
					ct = ct->parent();
				}
			}

			// Create cluster structure
			// For each vertex in the indyBag we create the cluster path
			// to the bag root if necessary.

			// Now build the cluster structure top down
			// Lists of root are never both empty
			List<cluster> queue;
			queue.pushBack(ca.indyBagRoot(i));
			cCopy[queue.front()] = bagCG.rootCluster();
			while (!queue.empty())
			{
				cluster c = queue.popFrontRet();

				//vertices are assigned to root by construction
				if (cCopy[c] != bagCG.rootCluster())
				{
					for(node u : cNodes[c])
						bagCG.reassignNode(u, cCopy[c]);
				}

				for(cluster ci : cChildren[c])
				{
					cCopy[ci] = bagCG.newCluster(cCopy[c]);
					queue.pushBack(ci);
				}

			}
#ifdef OGDF_DEBUG
			Logger::slout() << "Created clustered graph for indy bag with "<<theGraphs[i]->numberOfNodes()<< " nodes and "<<
					bagCG.numberOfClusters()<<" clusters\n";
			// Make sure the cluster structure is a rooted tree
			cluster t = bagCG.rootCluster();
			int ccnt = 0;
			List<cluster> cqueue;
			cqueue.pushBack(t);
			while (!cqueue.empty())
			{
				t = cqueue.popFrontRet();
				for(cluster c : t->children) {
					cqueue.pushBack(c);
				}
				ccnt++;
			}
			OGDF_ASSERT(ccnt == bagCG.numberOfClusters());
			string filename = string("IndySubcgraph") + to_string(i) + ".gml";

			ClusterGraphAttributes CGA(bagCG);
			GraphIO::writeGML(CGA, filename);
#endif
			//now the actual test, similar to the one below...
			if (theGraphs[i]->numberOfNodes() > 2) //could even start at 4...
			{
				makeParallelFreeUndirected(*theGraphs[i]); //Could do this while creating copy
				Logger::slout()<< "IndyBag of size n m c"<<theGraphs[i]->numberOfNodes()<< " "<<
						theGraphs[i]->numberOfEdges()<< " "<< bagCG.numberOfClusters()<<"\n";
				List<nodePair> ae;
				//Todo: Add an interface here that allows to transfer bag and
				// activity information to the master, otherwise we have to
				// compute this info twice.
				bool imresult = doTest(bagCG, ae);
	#ifdef OGDF_DEBUG
				Logger::slout() << "IndyBag number "<<i<<" is "<< (imresult ? "" : "non-") <<"c-planar\n";
				Logger::slout() << "Number of edges added for IndyBag: "<<ae.size()<<"\n";
	#endif
				result = result && imresult;
				if (!result) return result;
				for(const nodePair &np : ae) {
					addedEdges.emplaceBack(theGraphs[i]->original(np.v1), theGraphs[i]->original(np.v2));
				}
			}
#ifdef OGDF_DEBUG
			else
			{
				Logger::slout() << "IndyBag number "<<i<<" skipped due to size\n";
			}
#endif
		}//for indy bags

		for (int i = 0; i < numIndyBags; i++)
		{
			delete theGraphs[i];
		}
		delete [] theGraphs;

		// We test consistency by summing up the number of vertices.
	}
	else { //todo: can be joined again, as it is a special case wo cluster analysis
		//otherwise we just make a copy of the whole graph
		Graph G;
		ClusterArray<cluster> clusterCopy(CG);
		NodeArray<node> nodeCopy(CG.constGraph());
		EdgeArray<edge> edgeCopy(CG.constGraph());
		ClusterGraph C(CG,G,clusterCopy, nodeCopy, edgeCopy);
		makeParallelFreeUndirected(G);
		NodeArray<node> nodeOrig(G);
		for(node v : CG.constGraph().nodes)
		{
			nodeOrig[nodeCopy[v]] = v;
		}

		//Could use same list here for both graphs.
		addedEdges.clear();
		List<nodePair> ae;
		result = doTest(C,ae);
		//nodepairs are for the copy, store original nodes here
		for (const nodePair &np : ae) {
			addedEdges.emplaceBack(nodeOrig[np.v1], nodeOrig[np.v2]);
		}
	}

	return result;
}

bool ClusterPlanarity::doTest(const ClusterGraph &CG)
{
	List<nodePair> addEdges;
	return doTest(CG, addEdges);
}

bool ClusterPlanarity::doTest(const ClusterGraph &G,
			List<nodePair> &addedEdges)
{
//	if (m_solmeth==sm_new) return doFastTest(G,addedEdges);
	// We could take care of multiedges, but as long this is
	// not done, we do not allow this.
	OGDF_ASSERT(isParallelFreeUndirected(G)); // Graph has to be simple
#ifdef OGDF_DEBUG
	cout << "Creating new Masterproblem for clustergraph with "<<G.constGraph().numberOfNodes()<<" nodes\n";
#endif
	CP_MasterBase* cplanMaster;
	cplanMaster = new CPlanarityMaster(G,
									 m_heuristicLevel,
									 m_heuristicRuns,
									 m_heuristicOEdgeBound,
									 m_heuristicNPermLists,
									 m_kuratowskiIterations,
									 m_subdivisions,
									 m_kSupportGraphs,
									 m_kuratowskiHigh,
									 m_kuratowskiLow,
									 m_perturbation);
	if (m_solmeth == sm_new)
		static_cast<CPlanarityMaster*>(cplanMaster)->setSearchSpaceShrinking(true);
	else
		static_cast<CPlanarityMaster*>(cplanMaster)->setSearchSpaceShrinking(false);
		//new CPlanarMaster(G,m_heuristicLevel,m_heuristicRuns,m_heuristicOEdgeBound,m_heuristicNPermLists,m_kuratowskiIterations,
		//m_subdivisions,m_kSupportGraphs,m_kuratowskiHigh, m_kuratowskiLow,m_perturbation,m_branchingGap,m_time, m_pricing,
		//m_numAddVariables,m_strongConstraintViolation,m_strongVariableViolation,m_ol);
	cplanMaster->setTimeLimit(m_time.c_str());
	cplanMaster->setPortaFile(m_portaOutput);
	cplanMaster->useDefaultCutPool() = m_defaultCutPool;
#ifdef OGDF_DEBUG
	cout << "Starting Optimization\n";
#endif
	Master::STATUS abastatus;
	try {
		abastatus = cplanMaster->optimize();
	}
	catch (...)
	{
		#ifdef OGDF_DEBUG
		cerr << "ABACUS Optimization failed...\n";
		#endif
	}

	m_optStatus    = abastatus;
	m_totalTime    = getDoubleTime(*cplanMaster->totalTime());
	m_heurTime     = getDoubleTime(*cplanMaster->improveTime());
	m_sepTime      = getDoubleTime(*cplanMaster->separationTime());
	m_lpTime       = getDoubleTime(*cplanMaster->lpTime());
	m_lpSolverTime = getDoubleTime(*cplanMaster->lpSolverTime());
	m_totalWTime   = getDoubleTime(*cplanMaster->totalCowTime());
	m_numKCons     = cplanMaster->addedKConstraints();
	m_numCCons     = cplanMaster->addedCConstraints();
	m_numLPs       = cplanMaster->nLp();
	m_numBCs       = cplanMaster->nSub();
	m_numSubSelected = cplanMaster->nSubSelected();
	m_numVars      = cplanMaster->nMaxVars()-cplanMaster->getNumInactiveVars();
#ifdef OGDF_DEBUG
	m_solByHeuristic = cplanMaster->m_solByHeuristic;
#endif
#ifdef OGDF_DEBUG
	if(cplanMaster->pricing())
		Logger::slout() << "Pricing was ON\n";
	Logger::slout()<<"ABACUS returned with status '"<< Master::STATUS_[abastatus] <<"'\n"<<flush;
#endif

	cplanMaster->getConnectionOptimalSolutionEdges(addedEdges);
	//int addE = addedEdges.size();

#ifdef OGDF_DEBUG
	cout<<"-Number of added edges "<< addedEdges.size()<<"\n";
#endif

	if (m_portaOutput)
	{
		writeFeasible(getPortaFileName(), *cplanMaster, abastatus);
	}

	CP_MasterBase::solutionstate status = cplanMaster->m_solState;

	delete cplanMaster;

	switch (status) {
	case CP_MasterBase::ss_cp: return true; break;
	//Todo: catch and publish errors here
	case CP_MasterBase::ss_ncp: return false; break;
	default: cerr<<"** Undefined optimization result for c-planarity computation **\n"; break;
	}//switch

	return false; //Todo: Throw error here if eg outofmemory etc
}//doTest


//returns list of all clusters in subtree at c in bottom up order
void ClusterPlanarity::getBottomUpClusterList(const cluster c, List< cluster > & theList)
{
	for(cluster cc : c->children) {
		getBottomUpClusterList(cc, theList);
	}
	theList.pushBack(c);
}

//outputs the set of feasible solutions
void ClusterPlanarity::writeFeasible(const char *filename,
	CP_MasterBase &master,
	Master::STATUS &status)
{
	const ClusterGraph& CG = *(master.getClusterGraph());
	const Graph& G = CG.constGraph();
	//first compute the nodepairs that are potential candidates to connect
	//chunks in a cluster
	//potential connection edges
	NodeArray< NodeArray<bool> > potConn(G);
	for(node v : G.nodes)
	{
		potConn[v].init(G, false);
	}
	//we perform a bottom up cluster tree traversal
	List< cluster > clist;
	getBottomUpClusterList(CG.rootCluster(), clist);
	//could use postordertraversal instead

	List< nodePair > connPairs; //holds all connection node pairs
	//counts the number of potential connectivity edges
	//int potCount = 0; //equal to number of true values in potConn

	//we run through the clusters and check connected components
	//we consider all possible edges connecting CCs in a cluster,
	//even if they may be connected by edges in a child cluster
	//(to get the set of all feasible solutions)

	for(cluster c : clist)
	{
		//we compute the subgraph induced by vertices in c
		GraphCopy gcopy;
		gcopy.createEmpty(G);
		List<node> clusterNodes;
		//would be more efficient if we would just merge the childrens' vertices
		//and add c's
		c->getClusterNodes(clusterNodes);
		NodeArray<bool> activeNodes(G, false); //true for all cluster nodes
		EdgeArray<edge> copyEdge(G); //holds the edge copy

		for(node v : clusterNodes)
			activeNodes[v] = true;

		gcopy.initByActiveNodes(clusterNodes, activeNodes, copyEdge);
		//gcopy now represents the cluster induced subgraph

		//we compute the connected components and store all nodepairs
		//that connect two of them
		NodeArray<int> component(gcopy);
		connectedComponents(gcopy, component);
		//now we run over all vertices and compare the component
		//number of adjacent vertices. If they differ, we found a
		//potential connection edge. We do not care if we find them twice.
		for(node v : gcopy.nodes)
		{
			for(node w : gcopy.nodes)
			{
				if (component[v] != component[w])
				{
					cout <<"Indizes: "<<v->index()<<":"<<w->index()<<"\n";
					node vg = gcopy.original(v);
					node wg = gcopy.original(w);
					bool newConn = !((vg->index() < wg->index()) ? potConn[vg][wg] : potConn[wg][vg]);
					if (newConn)
					{
						nodePair np; np.v1 = vg; np.v2 = wg;
						connPairs.pushBack(np);
						if (vg->index() < wg->index())
							potConn[vg][wg] = true;
						else
							potConn[wg][vg] = true;
					}
				}
			}//nodes
		}//nodes
	}

	cout << "Number of potential connection edges: "<< connPairs.size()<<"\n";

	//we run through our candidates and save them in an array
	//that can be used for dynamic graph updates
	int i = 0;
	connStruct *cons = new connStruct[connPairs.size()];
	for(const nodePair &np : connPairs)
	{
		connStruct cs;
		cs.connected = false;
		cs.v1 = np.v1;
		cs.v2 = np.v2;
		cs.e  = nullptr;

		cons[i] = cs;
		i++;
	}

	//-------------------------------------------------------------------------
	// WARNING: this is extremely slow for graphs with a large number of cluster
	// chunks now we test all possible connection edge combinations for c-planarity
	Graph G2;

	NodeArray<node> origNodes(CG.constGraph());
	ClusterArray<cluster> origCluster(CG);
	EdgeArray<edge> origEdges(CG.constGraph());
	ClusterGraph testCopy(CG, G2, origCluster, origNodes, origEdges);

	ofstream os(filename);

	// Output dimension of the LP (number of variables)
	os << "DIM = " << connPairs.size() << "\n";
	os << "COMMENT\n";

	switch (status) {
		case Master::Optimal:
			os << "Optimal \n\n"; break;
		case Master::Error:
			os << "Error \n\n"; break;
		default:
			os << "unknown \n\n";
	}

	for (i = 0; i < connPairs.size(); i++)
	{
		os << "Var " << i << ": " << origNodes[cons[i].v1]->index() << "->" << origNodes[cons[i].v2] << "\n";
	}

	os << "CONV_SECTION\n";

#ifdef writefeasiblegraphs
	int j = 0; //debug
#endif

	if (connPairs.size() > 0)
	while (true)
	{
		//we create the next test configuration by incrementing the edge selection array
		//we create the corresponding graph dynamically on the fly
		i = 0;
		while ( (i < connPairs.size()) && (cons[i].connected == true) )
		{
			cons[i].connected = false;
			OGDF_ASSERT(cons[i].e != 0);
			G2.delEdge(cons[i].e);
			i++;
		}//while
		if (i >= connPairs.size()) break;
		//cout<<"v1graph: "<<&(*(cons[i].v1->graphOf()))<<"\n";
		//cout<<"origNodesgraph: "<<&(*(origNodes.graphOf()))<<"\n";
		cons[i].connected = true; //i.e., (false) will never be a feasible solution
		cons[i].e = G2.newEdge(origNodes[cons[i].v1], origNodes[cons[i].v2]);


		//and test it for c-planarity
		CconnectClusterPlanar CCCP;
		bool cplanar = CCCP.call(testCopy);

		//c-planar graphs define a feasible solution
		if (cplanar)
		{
#ifdef OGDF_DEBUG
			cout << "Feasible solution found\n";
#endif
			for (int j = 0; j < connPairs.size(); j++)
			{
				char ch = (cons[j].connected ? '1' : '0');
				cout << ch;
				os << ch << " ";
			}
			cout << "\n";
			os << "\n";
#ifdef writefeasiblegraphs
			string fn = "cGraph";
			fn += to_string(j++) + ".gml";
			GraphIO::writeGML(testCopy, fn);
#endif
		}
	}//while counting

	delete[] cons;

	os << "\nEND" <<"\n";
	os.close();

	//return;

	os.open(getIeqFileName());
	os << "DIM = " << m_numVars << "\n";
	// Output the status as a comment
	os << "COMMENT\n";

	switch (status) {
		case Master::Optimal:
			os << "Optimal \n\n"; break;
		case Master::Error:
			os << "Error \n\n"; break;
		default:
			os << "unknown \n\n";
	}

	// In case 0 is not a valid solution, some PORTA functions need
	//a valid solution in the ieq file
	os << "VALID\n";

	os << "\nLOWER_BOUNDS\n";

	for (i = 0; i < m_numVars; i++) os << "0 ";
	os << "\n";

	os << "\nHIGHER_BOUNDS\n";
	for (i = 0; i < m_numVars; i++) os << "1 ";
	os << "\n";

	os << "\nINEQUALITIES_SECTION\n";
	//we first read the standard constraint that are written
	//into a text file by the optimization master
	ifstream isf(master.getStdConstraintsFileName());
	if (!isf)
	{
		cerr << "Could not open optimization master's standard constraint file\n";
		os << "#No standard constraints read\n";
	}
	else
	{
		char* fileLine = new char[maxConLength()];
		while (isf.getline(fileLine, maxConLength()))
		{
			//skip commment lines
			if (fileLine[0] == '#') continue;
			int count = 1;
			std::istringstream iss(fileLine);
			char d;
			bool rhs = false;
			while (iss >> d)
			{
				if ( rhs || ( (d == '<') || (d == '>') || (d == '=') ) )
				{
					os << d;
					rhs = true;
				}
				else
				{
					if (d != '0')
					{
						os <<"+"<< d <<"x"<<count;
					}
					count++;
				}
			}//while chars
			os <<"\n";
		}
		delete[] fileLine;
	}//ifstream
	//now we read the cut pools from the master
	if (master.useDefaultCutPool())
	{
		os << "#No cut constraints read from master\n";
		//StandardPool<Constraint, Variable> *connCon = master.cutPool();
	}
	else
	{
		StandardPool<Constraint, Variable> *connCon = master.getCutConnPool();
		StandardPool<Constraint, Variable> *kuraCon = master.getCutKuraPool();
		StandardPool<Variable, Constraint> *stdVar = master.varPool();
		OGDF_ASSERT(connCon != 0);
		OGDF_ASSERT(kuraCon != 0);
		cout << connCon->number() << " Constraints im MasterConnpool \n";
		cout << kuraCon->number() << " Constraints im MasterKurapool \n";
		cout << connCon->size() << " Größe ConnPool"<<"\n";
		outputCons(os, connCon, stdVar);
		outputCons(os, kuraCon, stdVar);
	}//else
	os << "\nEND" <<"\n";
	os.close();
	cout << "Cutting is set: "<<master.cutting()<<"\n";
	//cout <<"Bounds for the variables:\n";
	//Sub &theSub = *(master.firstSub());
	//for ( i = 0; i < theSub.nVar(); i++)
	//{
	//	cout << i << ": " << theSub.lBound(i) << " - " << theSub.uBound(i) << "\n";
	//}
	/*// OLD CRAP
	cout << "Constraints: \n";
	StandardPool< Constraint, Variable > *spool = master.conPool();
	StandardPool< Constraint, Variable > *cpool = master.cutPool();

	cout << spool->size() << " Constraints im Masterpool \n";
	cout << cpool->size() << " Constraints im Mastercutpool \n";

	cout << "ConPool Constraints \n";
	for ( i = 0; i < spool->size(); i++)
	{
		PoolSlot< Constraint, Variable > * sloty = spool->slot(i);
		Constraint *mycon = sloty->conVar();
		switch (mycon->sense()->sense())
		{
			case CSense::Less: cout << "<" << "\n"; break;
			case CSense::Greater: cout << ">" << "\n"; break;
			case CSense::Equal: cout << "=" << "\n"; break;
			default: cout << "Inequality sense doesn't make any sense \n"; break;
		}//switch
	}
	cout << "CutPool Constraints \n";
	for ( i = 0; i < cpool->size(); i++)
	{
		PoolSlot< Constraint, Variable > * sloty = cpool->slot(i);
		Constraint *mycon = sloty->conVar();
		switch (mycon->sense()->sense())
		{
			case CSense::Less: cout << "<" << "\n"; break;
			case CSense::Greater: cout << ">" << "\n"; break;
			case CSense::Equal: cout << "=" << "\n"; break;
			default: cout << "Inequality sense doesn't make any sense \n"; break;
		}//switch
	}
	*/
	/*
	for ( i = 0; i < theSub.nCon(); i++)
	{
		Constraint &theCon = *(theSub.constraint(i));

		for ( i = 0; i < theSub.nVar(); i++)
		{
			double c = theCon.coeff(theSub.variable(i));
			if (c != 0)
				cout << c;
			else cout << "  ";
		}
		switch (theCon.sense()->sense())
		{
			case CSense::Less: cout << "<" << "\n"; break;
			case CSense::Greater: cout << ">" << "\n"; break;
			case CSense::Equal: cout << "=" << "\n"; break;
			default: cout << "doesn't make any sense \n"; break;
		}//switch
		float fl;
		while(!(std::cin >> fl))
		{
			std::cin.clear();
			std::cin.ignore(numeric_limits<streamsize>::max(),'\n');
		}
	}*/
}//writeportaieq

void ClusterPlanarity::outputCons(ofstream &os,
	StandardPool<Constraint, Variable> *connCon,
	StandardPool<Variable, Constraint> *stdVar)
{
	int i;
	for ( i = 0; i < connCon->number(); i++)
		{
			PoolSlot< Constraint, Variable > * sloty = connCon->slot(i);
			Constraint *mycon = sloty->conVar();
			OGDF_ASSERT(mycon != 0);
			int count;
			for (count = 0; count < stdVar->size(); count++)
			{
				PoolSlot< Variable, Constraint > * slotv = stdVar->slot(count);
				Variable *myvar = slotv->conVar();
				double d = mycon->coeff(myvar);
				if (d != 0.0) //precision!
				{
					os <<"+"<< d <<"x"<<count+1;
				}
			}//for
			switch (mycon->sense()->sense())
			{
				case CSense::Less: os << " <= "; break;
				case CSense::Greater: os << " >= "; break;
				case CSense::Equal: os << " = "; break;
				default: os << "Inequality sense doesn't make any sense \n";
						 cerr << "Inequality sense unknown \n";
						break;
			}//switch
			os << mycon->rhs();
			os << "\n";
		}
}

} //end namespace ogdf

#endif //USE_ABACUS
