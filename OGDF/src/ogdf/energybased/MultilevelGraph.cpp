/** \file
 * \brief MLG is the main data structure for ModularMultilevelMixer
 *
 * \author Gereon Bartel
 *
 * \par License:
 * This file is part of the Open Graph Drawing Framework (OGDF).
 *
 * \par
 * Copyright (C)<br>
 * See README.txt in the root directory of the OGDF installation for details.
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

#include <ogdf/internal/energybased/MultilevelGraph.h>
#include <ogdf/basic/simple_graph_alg.h>
#include <ogdf/fileformats/GraphIO.h>



namespace ogdf {


MultilevelGraph::~MultilevelGraph()
{
	// delete all Nodemerges!
	while(!m_changes.empty()) {
		delete m_changes.back();
		m_changes.pop_back();
	}

	delete m_GA;
	m_reverseNodeIndex.clear();
	//only delete the Graph if it was created!
	if (m_createdGraph)	{
		delete m_G;
	}

	m_reverseEdgeIndex.clear();
}

//initialize internal structures such as the GraphAttributes that store the layout
void MultilevelGraph::initInternal()
{
	OGDF_ASSERT(m_G != 0);
	m_GA = new GraphAttributes(*m_G);
}

MultilevelGraph::MultilevelGraph()
:m_createdGraph(true)
{
	m_G = new Graph();
	if(m_G == nullptr) {
		OGDF_THROW(InsufficientMemoryException);
	}

	//replaces layout info stuff below
	initInternal();

	m_nodeAssociations.init(*m_G, 0);
	m_edgeAssociations.init(*m_G, 0);
	m_radius.init(*m_G, 1.0);
	m_weight.init(*m_G, 1.0);

	initReverseIndizes();
}


MultilevelGraph::MultilevelGraph(GraphAttributes &GA)
:m_createdGraph(true)
{
	m_G = new Graph();
	if(m_G == nullptr) {
		OGDF_THROW(InsufficientMemoryException);
	}

	//replaces layout info stuff below
	initInternal();

	m_nodeAssociations.init(*m_G);
	m_edgeAssociations.init(*m_G);
	m_radius.init(*m_G);
	m_weight.init(*m_G);
	copyFromGraph(GA.constGraph(), m_nodeAssociations, m_edgeAssociations);
	prepareGraphAttributes(GA);
	importAttributes(GA);

	//initReverseIndizes();
}


MultilevelGraph::MultilevelGraph(Graph &G)
:m_createdGraph(false), m_G(nullptr)
{
	m_G = &G;

	//replaces layout info stuff below
	initInternal();

	m_nodeAssociations.init(*m_G, 0);
	m_edgeAssociations.init(*m_G, 0);
	m_radius.init(*m_G, 1.0);
	m_weight.init(*m_G, 1.0);

	initReverseIndizes();
}


MultilevelGraph::MultilevelGraph(GraphAttributes &GA, Graph &G)
:m_createdGraph(false), m_G(nullptr)
{
	m_G = &G;
	m_nodeAssociations.init(*m_G, 0);
	m_edgeAssociations.init(*m_G, 0);
	m_radius.init(*m_G);
	m_weight.init(*m_G);

	initInternal();

	prepareGraphAttributes(GA);
	importAttributes(GA);

	initReverseIndizes();
}


MultilevelGraph::MultilevelGraph(istream &is)
:m_createdGraph(true)
{
	m_G = new Graph();
	if(m_G == nullptr) {
		OGDF_THROW(InsufficientMemoryException);
	}
	m_nodeAssociations.init(*m_G);
	m_edgeAssociations.init(*m_G);
	m_radius.init(*m_G);
	m_weight.init(*m_G);

	initInternal();
	//GraphAttributes tempGA(*m_G);
	GraphIO::readGML(*m_GA, *m_G, is);
	prepareGraphAttributes(*m_GA);
	importAttributesSimple(*m_GA);

	initReverseIndizes();
}


MultilevelGraph::MultilevelGraph(const char *filename) : m_createdGraph(true)
{
	m_G = new Graph();
	if(m_G == nullptr) {
		OGDF_THROW(InsufficientMemoryException);
	}
	m_nodeAssociations.init(*m_G);
	m_edgeAssociations.init(*m_G);
	m_radius.init(*m_G);
	m_weight.init(*m_G);

	initInternal();
	//GraphAttributes tempGA(*m_G);
	GraphIO::readGML(*m_GA, *m_G, filename);
	prepareGraphAttributes(*m_GA);
	importAttributesSimple(*m_GA);

	initReverseIndizes();
}


void MultilevelGraph::prepareGraphAttributes(GraphAttributes &GA) const
{
	long additionalAttributes = 0;
	if (!(GA.attributes() & GraphAttributes::edgeDoubleWeight)) {
		additionalAttributes |= GraphAttributes::edgeDoubleWeight;
	}
	if (!(GA.attributes() & GraphAttributes::nodeWeight)) {
		additionalAttributes |= GraphAttributes::nodeWeight;
	}
	GA.initAttributes(additionalAttributes);
}


void MultilevelGraph::copyFromGraph(const Graph &G, NodeArray<int> & /*nodeAssociations*/, EdgeArray<int> & /* edgeAssociations */)
{
	NodeArray<node> tempAssociations(G);

	for(node v : G.nodes) {
		node v_new = m_G->newNode();
		m_nodeAssociations[v_new] = v->index();
		tempAssociations[v] = v_new;
	}

	for(edge e : G.edges) {
		edge e_new = m_G->newEdge(tempAssociations[e->source()], tempAssociations[e->target()]);
		m_edgeAssociations[e_new] = e->index();
	}

	initReverseIndizes();
}


int MultilevelGraph::getLevel()
{
	if (m_changes.size() == 0)
	{
		return 0;
	}
	else
	{
		return m_changes.back()->m_level;
	}
}


// assumes, that the Graphs of MultilevelGraph and GA are the same, not copies!
void MultilevelGraph::exportAttributesSimple(GraphAttributes &GA) const
{
	OGDF_ASSERT(&(GA.constGraph()) == m_G);

	prepareGraphAttributes(GA);

	for(node v : m_G->nodes) {
		GA.x(v) =  m_GA->x(v);
		GA.y(v) =  m_GA->y(v);
		//TODO: Check what this w,h computation does
		double w = GA.width(v);
		double h = GA.height(v);
		if(w > 0 || h > 0) {
			double factor =  m_radius[v] / sqrt(w*w + h*h) * 2.0f;
			w *= factor;
			h *= factor;
		} else {
			w = h = m_radius[v] * sqrt(2.0f);
		}
		GA.width(v) = w;
		GA.height(v) = h;
		GA.weight(v) = m_reverseNodeMergeWeight[v->index()];
	}

	for(edge e : m_G->edges) {
		GA.doubleWeight(e) = m_weight[e];
	}
}


void MultilevelGraph::exportAttributes(GraphAttributes &GA) const
{
	OGDF_ASSERT(GA.constGraph().numberOfNodes() == m_G->numberOfNodes());
	OGDF_ASSERT(GA.constGraph().numberOfEdges() == m_G->numberOfEdges());

	prepareGraphAttributes(GA);

	std::vector<node> tempNodeAssociations;
	const Graph &cG = GA.constGraph();
	tempNodeAssociations.resize(cG.maxNodeIndex()+1, nullptr);

	for(node v : cG.nodes) {
		tempNodeAssociations[v->index()] = v;
	}

	for(node v : m_G->nodes) {
		GA.x(tempNodeAssociations[m_nodeAssociations[v]]) =  m_GA->x(v);
		GA.y(tempNodeAssociations[m_nodeAssociations[v]]) =  m_GA->y(v);
		double w = GA.width(tempNodeAssociations[m_nodeAssociations[v]]);
		double h = GA.height(tempNodeAssociations[m_nodeAssociations[v]]);
		if(w > 0 || h > 0) {
			double factor =  m_radius[v] / sqrt(w*w + h*h) * 2.0f;
			w *= factor;
			h *= factor;
		} else {
			w = h = m_radius[v] * sqrt(2.0f);
		}
		GA.width(tempNodeAssociations[m_nodeAssociations[v]]) = w;
		GA.height(tempNodeAssociations[m_nodeAssociations[v]]) = h;
		GA.weight(tempNodeAssociations[m_nodeAssociations[v]]) = m_reverseNodeMergeWeight[v->index()];
	}

	std::vector<edge> tempEdgeAssociations;
	tempEdgeAssociations.resize(cG.maxEdgeIndex()+1, nullptr);
	for(edge e :cG.edges) {
		tempEdgeAssociations[e->index()] = e;
	}

	for(edge e : m_G->edges) {
		GA.doubleWeight(tempEdgeAssociations[m_edgeAssociations[e]]) = m_weight[e];
	}
}


void MultilevelGraph::importAttributesSimple(const GraphAttributes &GA)
{
	OGDF_ASSERT(&(GA.constGraph()) == m_G);

	m_avgRadius = 0.0;

	for(node v : m_G->nodes) {
		double w = GA.width(v);
		double h = GA.height(v);
		if(w > 0 || h > 0) {
			m_radius[v] = sqrt(w*w + h*h) / 2.0f;
		} else {
			m_radius[v] = 1.0f;
		}
		m_avgRadius += m_radius[v];
		m_GA->x(v) = GA.x(v);
		m_GA->y(v) = GA.y(v);
		m_GA->width(v) = GA.width(v);
		m_GA->height(v) = GA.height(v);
	}
	m_avgRadius /= m_G->numberOfNodes();

	for(edge e : m_G->edges) {
		m_weight[e] = GA.doubleWeight(e);
	}
}


void MultilevelGraph::importAttributes(const GraphAttributes &GA)
{
	OGDF_ASSERT(GA.constGraph().numberOfNodes() == m_G->numberOfNodes());
	OGDF_ASSERT(GA.constGraph().numberOfEdges() == m_G->numberOfEdges());

	m_avgRadius = 0.0;

	std::vector<node> tempNodeAssociations;
	const Graph &cG = GA.constGraph();
	tempNodeAssociations.resize(cG.maxNodeIndex()+1, nullptr);
	for(node v : cG.nodes) {
		tempNodeAssociations[v->index()] = v;
	}

	for(node v : m_G->nodes) {

		double w = GA.width(tempNodeAssociations[m_nodeAssociations[v]]);
		double h = GA.height(tempNodeAssociations[m_nodeAssociations[v]]);
		if(w > 0 || h > 0) {
			m_radius[v] = sqrt(w*w + h*h) / 2.0f;
		} else {
			m_radius[v] = 1.0f;
		}

		m_avgRadius += m_radius[v];

		m_GA->x(v) = GA.x(tempNodeAssociations[m_nodeAssociations[v]]);
		m_GA->y(v) = GA.y(tempNodeAssociations[m_nodeAssociations[v]]);
		m_GA->width(v) = GA.width(tempNodeAssociations[m_nodeAssociations[v]]);
		m_GA->height(v) = GA.height(tempNodeAssociations[m_nodeAssociations[v]]);
	}

	m_avgRadius /= m_G->numberOfNodes();

	std::vector<edge> tempEdgeAssociations;
	tempEdgeAssociations.resize(cG.maxEdgeIndex()+1, nullptr);
	for(edge e : cG.edges) {
		tempEdgeAssociations[e->index()] = e;
	}

	for(edge e : m_G->edges) {
		m_weight[e] = GA.doubleWeight(tempEdgeAssociations[m_edgeAssociations[e]]);
	}
}


void MultilevelGraph::reInsertGraph(MultilevelGraph &MLG)
{
	std::map<node, node> tempNodeAssociations;

	for(node v : MLG.m_G->nodes) {
		MLG.copyNodeTo(v, *this, tempNodeAssociations, false, MLG.m_nodeAssociations[v]);
	}

	for(edge e : MLG.m_G->edges) {
		MLG.copyEdgeTo(e, *this, tempNodeAssociations, false, MLG.m_edgeAssociations[e]);
	}

	initReverseIndizes();
}


void MultilevelGraph::reInsertAll(std::vector<MultilevelGraph *> &components)
{
	for(MultilevelGraph *g : components)
	{
		reInsertGraph(*g);
	}
}


// keeps Changes
// keeps Node and Edge Associations
// deletes Nodes and Eges from Graph
// deletes Attributes
// deprecated, use componentsplitterlayout instead
std::vector<MultilevelGraph *> MultilevelGraph::splitIntoComponents()
{
	std::vector<MultilevelGraph *> components;

	NodeArray<int> componentNumbers(*m_G);
	int numComponents = connectedComponents(*m_G, componentNumbers);
	if (numComponents == 0) {
		return components;
	}

	std::vector< std::vector<node> > componentArray;
	componentArray.resize(numComponents);
	for(node v : m_G->nodes) {
		componentArray[componentNumbers[v]].push_back(v);
	}

	for (unsigned int componentNumber = 0; componentNumber < componentArray.size(); componentNumber++) {
		std::vector<node> componentSubArray = componentArray[componentNumber];
		MultilevelGraph * component = removeOneCC(componentSubArray);
		components.push_back(component);
	}

	OGDF_ASSERT(m_G->numberOfNodes() == 0);
	OGDF_ASSERT(m_G->numberOfEdges() == 0);

	m_radius.init(*m_G);
	m_weight.init(*m_G);

	return components;
}


void MultilevelGraph::copyNodeTo(node v, MultilevelGraph &MLG, std::map<node, node> &tempNodeAssociations, bool associate, int index)
{
	node v_new;
	if (index == -1) {
		v_new = MLG.m_G->newNode();
	} else {
		v_new = MLG.m_G->newNode(index);
	}

	tempNodeAssociations[v] = v_new;
	if(associate) {
		MLG.m_nodeAssociations[v_new] = v->index();
	}
	MLG.m_radius[v_new] = m_radius[v];
	MLG.x(v_new, x(v));
	MLG.y(v_new, y(v));
}


void MultilevelGraph::copyEdgeTo(edge e, MultilevelGraph &MLG, std::map<node, node> &tempNodeAssociations, bool associate, int index)
{
	node source = e->source();
	node target = e->target();
	edge e_new;
	if (index == -1) {
		e_new = MLG.m_G->newEdge(tempNodeAssociations[source], tempNodeAssociations[target]);
	} else {
		e_new = MLG.m_G->newEdge(tempNodeAssociations[source], tempNodeAssociations[target], index);
	}

	if(associate) {
		MLG.m_edgeAssociations[e_new] = e->index();
	}
	MLG.m_weight[e_new] = m_weight[e];
}


MultilevelGraph * MultilevelGraph::removeOneCC(std::vector<node> &componentSubArray)
{
	MultilevelGraph * MLGcomponent = new MultilevelGraph();
	std::map<node, node> tempNodeAssociations;

	// copy nodes
	for (node v : componentSubArray) {
		copyNodeTo(v, *MLGcomponent, tempNodeAssociations, true);
	}

	// move edges
	for (node v : componentSubArray) {
		edge e;
//		std::vector<edge> toDelete;
		forall_adj_edges(e, v) {
			if (e != nullptr && e->source() == v) {
				copyEdgeTo(e, *MLGcomponent, tempNodeAssociations, true);
//				toDelete.push_back(e);
			}
		}
/*		// Test if this is good for Performace.
		// makes Assert Edges == 0 fail!
		// Because of self loops!
		for(std::vector<edge>::iterator j = toDelete.begin(); j != toDelete.end(); j++) {
			m_G->delEdge(*j);
		}
*/
	}

	tempNodeAssociations.clear();

	// delete nodes
	for (node v : componentSubArray) {
		m_G->delNode(v);
	}

	MLGcomponent->initReverseIndizes();
	return MLGcomponent;
}


bool MultilevelGraph::postMerge(NodeMerge * NM, node merged)
{
	// merged has no more edges!
	int index = merged->index();
	if (merged->degree() == 0 && NM->m_changedNodes.size() > 0) {
		NM->m_mergedNode = index;
		NM->m_radius[index] = m_radius[index];
		m_changes.push_back(NM);
		m_G->delNode(merged);
		m_reverseNodeIndex[index] = nullptr;
		return true;
	} else {
		return false;
	}
}


bool MultilevelGraph::changeNode(NodeMerge * NM, node theNode, double newRadius, node merged)
{
	int index = theNode->index();
	//we assume that changeNode is called exactly onces when a node is merged
	//with its parent with parameter theNode being the parent and add 1 to
	//the parents merge weight
	m_reverseNodeMergeWeight[index] += m_reverseNodeMergeWeight[merged->index()];
	std::vector<int>::iterator pos = find(NM->m_changedNodes.begin(), NM->m_changedNodes.end(), index);

	if (pos == NM->m_changedNodes.end()) {
		NM->m_changedNodes.push_back(index);
		NM->m_radius[index] = m_radius[theNode];
	}
	m_radius[theNode] = newRadius;

	return true;
}


bool MultilevelGraph::changeEdge(NodeMerge * NM, edge theEdge, double newWeight, node newSource, node newTarget)
{
	int index = theEdge->index();
	std::vector<int>::iterator pos = find(NM->m_changedEdges.begin(), NM->m_changedEdges.end(), index);

	if (pos == NM->m_changedEdges.end()) {
		NM->m_changedEdges.push_back(index);
		NM->m_doubleWeight[index] = m_weight[theEdge];
		NM->m_source[index] = theEdge->source()->index();
		NM->m_target[index] = theEdge->target()->index();
	}
	m_G->delEdge(theEdge);
	m_reverseEdgeIndex[index] = m_G->newEdge(newSource, newTarget, index);
	m_weight[theEdge] = newWeight;

	return true;
}


bool MultilevelGraph::deleteEdge(NodeMerge * NM, edge theEdge)
{
	int index = theEdge->index();

	NM->m_deletedEdges.push_back(index);
	NM->m_doubleWeight[index] = m_weight[theEdge];
	NM->m_source[index] = theEdge->source()->index();
	NM->m_target[index] = theEdge->target()->index();

	m_G->delEdge(theEdge);
	m_reverseEdgeIndex[index] = nullptr;

	return true;
}


std::vector<edge> MultilevelGraph::moveEdgesToParent(NodeMerge * NM, node theNode, node parent, bool deleteDoubleEdges, int adjustEdgeLengths)
{
	OGDF_ASSERT(theNode != parent);

	std::vector<edge> doubleEdges;
	std::vector<edge> adjEdges;
	edge e;
	forall_adj_edges(e, theNode) {
		adjEdges.push_back(e);
	}

	double nodeToParentLen = 0.0;
	for (edge e : adjEdges)
	{
		node newSource = e->source();
		node newTarget = e->target();
		if ((newSource == theNode && newTarget == parent)
		|| (newSource == parent && newTarget == theNode)){
			nodeToParentLen = m_weight[e];
			break;
		}
	}

	for (edge e : adjEdges)
	{
		node newSource = e->source();
		node newTarget = e->target();

		if (newSource == theNode) {
			newSource = parent;
		}
		if (newTarget == theNode) {
			newTarget = parent;
		}

		bool exists = false;
		edge twinEdge = nullptr;
		for(adjEntry adj : parent->adjEdges) {
			if (adj->twinNode() != parent && (adj->twinNode() == newSource || adj->twinNode() == newTarget)) {
				exists = true;
				twinEdge = adj->theEdge();
				double extraLength = 0.0;
				if(adjustEdgeLengths != 0) {
					extraLength = m_weight[twinEdge] + adjustEdgeLengths * nodeToParentLen;
				}
				changeEdge(NM, twinEdge, (m_weight[twinEdge] + m_weight[e] + extraLength) * 0.5f, twinEdge->source(), twinEdge->target());
				break;
			}
		}

		// has this edge already
		if (exists || newSource == newTarget) {
			doubleEdges.push_back(e);
		} else {
			changeEdge(NM, e, m_weight[e], newSource, newTarget);
		}
	}

	if (deleteDoubleEdges) {
		while (!doubleEdges.empty()) {
			deleteEdge(NM, doubleEdges.back());
			doubleEdges.pop_back();
		}
	}

	OGDF_ASSERT(theNode->degree() == (int)doubleEdges.size());

	// not deleted edges that are adjacent to theNode are returned.
	return doubleEdges;
}


NodeMerge * MultilevelGraph::getLastMerge()
{
	return m_changes.back();
}


node MultilevelGraph::undoLastMerge()
{
	if (m_changes.empty()) {
		return nullptr;
	}
	NodeMerge * merge = m_changes.back();
	m_changes.pop_back();

	// reinsert merged node
	int index = merge->m_mergedNode;
	node merged = m_G->newNode(index);
	m_reverseNodeIndex[index] = merged;
	m_radius[merge->m_mergedNode] = merge->m_radius[index];

	// add deleted edges
	for (int index : merge->m_deletedEdges) {
		m_reverseEdgeIndex[index] = m_G->newEdge(m_reverseNodeIndex[merge->m_source[index]], m_reverseNodeIndex[merge->m_target[index]], index);
		m_weight[index] = merge->m_doubleWeight[index];
	}

	// undo edge changes
	for (int index : merge->m_changedEdges) {
		m_G->delEdge(m_reverseEdgeIndex[index]);
		m_reverseEdgeIndex[index] = m_G->newEdge(m_reverseNodeIndex[merge->m_source[index]], m_reverseNodeIndex[merge->m_target[index]], index);
		m_weight[index] = merge->m_doubleWeight[index];
	}

	// undo node changes
	for (int index : merge->m_changedNodes) {
		m_radius[index] = merge->m_radius[index];
		m_reverseNodeMergeWeight[index] -= m_reverseNodeMergeWeight[merged->index()];
	}

	delete merge;
	return merged;
}


edge MultilevelGraph::getEdge(unsigned int index)
{
	if (index >= m_reverseEdgeIndex.size()) {
		return nullptr;
	}
	return m_reverseEdgeIndex[index];
}


node MultilevelGraph::getNode(unsigned int index)
{
	if (index >= m_reverseNodeIndex.size()) {
		return nullptr;
	}
	return m_reverseNodeIndex[index];
}


void MultilevelGraph::initReverseIndizes()
{
	if (m_G->numberOfNodes() > 0) {
		m_reverseNodeIndex.resize(m_G->maxNodeIndex()+1, nullptr);
	}
	if (m_G->numberOfNodes() > 0) {
		m_reverseNodeMergeWeight.resize(m_G->maxNodeIndex()+1, 1);
	}
	if (m_G->numberOfEdges() > 0) {
		m_reverseEdgeIndex.resize(m_G->maxEdgeIndex()+1, nullptr);
	}
}


void MultilevelGraph::updateMergeWeights()
{
	for(node v : m_G->nodes) {
		m_reverseNodeMergeWeight[v->index()] = 1;
	}
}


void MultilevelGraph::updateReverseIndizes()
{
	if ((unsigned int)m_G->maxNodeIndex() >= m_reverseNodeIndex.size() || (unsigned int)m_G->maxEdgeIndex() >= m_reverseEdgeIndex.size()) {
		initReverseIndizes();
	}

	for(node v : m_G->nodes) {
		m_reverseNodeIndex[v->index()] = v;
	}

	for(edge e : m_G->edges) {
		m_reverseEdgeIndex[e->index()] = e;
	}
}


void MultilevelGraph::writeGML(ostream &os)
{
	GraphAttributes GA(*m_G);
	exportAttributesSimple(GA);

	GraphIO::writeGML(GA, os);
}


void MultilevelGraph::writeGML(const char *fileName)
{
	ofstream os(fileName);
	writeGML(os);
}


void MultilevelGraph::moveToZero()
{
	// move Graph to zero
	double avg_x = 0.0;
	double avg_y = 0.0;
	for(node v : getGraph().nodes) {
		avg_x += x(v);
		avg_y += y(v);
	}
	avg_x /= getGraph().numberOfNodes();
	avg_y /= getGraph().numberOfNodes();
	for(node v : getGraph().nodes) {
		x(v, x(v) - avg_x);
		y(v, y(v) - avg_y);
	}
}


} // namespace ogdf
