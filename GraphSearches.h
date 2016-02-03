#ifndef CREEK_INCISOR_VISITOR
#define CREEK_INCISOR_VISITOR


#include <boost/graph/depth_first_search.hpp>
#include <boost/shared_ptr.hpp>
#include "Types.h"

int recursion_depth;

template <class Edge, class Graph>
class CycleTerminator : public boost::dfs_visitor<>
{
    
private:
    std::vector<Edge> & removal_list;

public:
    CycleTerminator(std::vector<Edge> & _removal_list) : removal_list(_removal_list)
    {
        
    }
    
    void back_edge(Edge e, const Graph& g) {
        removal_list.push_back(e);
    }
};

int
leafBranchVisit(Graph & channel_graph, EdgeDescriptor e, int distance, int cutoff, std::set<EdgeDescriptor> & edge_removal_list, std::set<VertexDescriptor> & vertex_removal_list)
{
    VertexDescriptor u = boost::target(e, channel_graph);
    ++distance;
    
    if (distance > cutoff) return (-1);
    
    if (channel_graph[u].type > 1)
    {
        edge_removal_list.insert(e);
        //vertex_removal_list.push_back(u);
        return (distance);
    }
    
    typedef boost::graph_traits<Graph>::out_edge_iterator OutIt;
    OutIt e_it, e_end;
    for (boost::tie(e_it, e_end) = boost::out_edges(u, channel_graph); e_it != e_end; e_it++)
    {
        if (*e_it != e)
        {
            int leaf_branch_dist = leafBranchVisit(channel_graph, *e_it, distance, cutoff, edge_removal_list, vertex_removal_list);
            if (leaf_branch_dist == -1) return (leaf_branch_dist);
            edge_removal_list.insert(e);
            vertex_removal_list.insert(u);
            return (leaf_branch_dist);
        }
    }
    return (-1);
}

void
leafBranchSearch(Graph& channel_graph, VertexDescriptor channel_vertex, int cutoff, std::set<EdgeDescriptor> & edge_removal_list, std::set<VertexDescriptor> & vertex_removal_list)
{
    typedef boost::graph_traits<Graph>::out_edge_iterator OutIt;
    OutIt e_it, e_end;
    int count = 0;
    int distance = 0;
    for (boost::tie(e_it, e_end) = boost::out_edges(channel_vertex, channel_graph); e_it != e_end; e_it++)
    {
        //There should only be one....
        distance = leafBranchVisit(channel_graph, *e_it, distance, cutoff, edge_removal_list, vertex_removal_list);
        ++count;
    }
    if (count > 1)
    {
        std::cerr << "Error: should only be one out edge of a terminal node" << std::endl;
    }
    if (distance != -1) vertex_removal_list.insert(channel_vertex);
}

//Made recursive version instead....
//void
//leafBranchRemoval(Graph& channel_graph, VertexDescriptor terminal_node, unsigned int trim_level, std::vector<EdgeDescriptor> & edge_removal_list, std::vector<VertexDescriptor> & vertex_removal_list)
//{
//    unsigned int branch_count = 0;
//    std::vector<EdgeDescriptor> loc_edge_removal_list;
//    std::vector<VertexDescriptor> loc_vertex_removal_list;
//    bool do_remove = false;
//    typedef boost::graph_traits<Graph>::out_edge_iterator OutIt;
//    OutIt e_it, e_end;
//    VertexDescriptor u = terminal_node;
//    EdgeDescriptor previousEdge;
//    EdgeDescriptor currentEdge;
//    while (branch_count < trim_level)
//    {
//        loc_vertex_removal_list.push_back(u);
//        for (boost::tie(e_it, e_end) = boost::out_edges(u, channel_graph); e_it != e_end; ++e_it)
//        {
//            if
//        };
//        
//        loc_edge_removal_list.push_back(*e_it);
//        VertexDescriptor u1 = boost::target(*e_it, channel_graph);
//        VertexDescriptor u2 = boost::source(*e_it, channel_graph);
//        previousEdge = *e_it;
//        if (u1 != u)
//        {
//            u = u1;
//        }
//        else if (u2 != u)
//        {
//            u = u2;
//        }
//        ++branch_count;
//        
//        if (channel_graph[u].type > 1)
//        {
//            do_remove = true;
//            break;
//        }
//    }
//    if (do_remove)
//    {
//        BOOST_FOREACH(EdgeDescriptor e, loc_edge_removal_list)
//        {
//            edge_removal_list.push_back(e);
//        }
//        BOOST_FOREACH(VertexDescriptor v, loc_vertex_removal_list)
//        {
//            vertex_removal_list.push_back(v);
//        }
//    }
//}



std::pair<VertexDescriptor, int>
controlSearchVisit(Graph & channel_graph, EdgeDescriptor e, int distance, VertexDescriptor start_vertex)
{
	++recursion_depth;
    VertexDescriptor u = boost::target(e, channel_graph);
    if (channel_graph[u].type > 1)
    {
        return (std::make_pair(u, distance));
    }
    else
    {
        if (channel_graph[u].up_cntrl_id == -1)
        {
            channel_graph[u].up_cntrl_id = channel_graph[start_vertex].node_id;
            channel_graph[u].up_cntrl_dist = distance;
        }
        else if (channel_graph[u].down_cntrl_id == -1)
        {
            channel_graph[u].down_cntrl_id = channel_graph[start_vertex].node_id;
            channel_graph[u].down_cntrl_dist = distance;
        }
        
        typedef boost::graph_traits<Graph>::out_edge_iterator OutIt;
        OutIt e_it, e_end;
        for (boost::tie(e_it, e_end) = boost::out_edges(u, channel_graph); e_it != e_end; e_it++)
        {
            if (*e_it != e)
            {
                ++distance;
                return (controlSearchVisit(channel_graph, *e_it, distance, start_vertex));
                break;
            }
        }
    }
    return (std::make_pair(u, distance));
}


void
controlSearch(Graph& channel_graph,  VertexDescriptor channel_vertex, Graph& controls_graph, VertexDescriptor control_vertex, VertexIDMap control_vertex_map)
{
    typedef boost::graph_traits<Graph>::out_edge_iterator OutIt;
    OutIt e_it, e_end;
    for (boost::tie(e_it, e_end) = boost::out_edges(channel_vertex, channel_graph); e_it != e_end; e_it++)
    {
        int distance = 1;
        VertexDescriptor channel_c;
		recursion_depth = 0;
        boost::tie(channel_c, distance) = controlSearchVisit(channel_graph, *e_it, distance, channel_vertex);
        VertexDescriptor cntrl_c = control_vertex_map[channel_graph[channel_c].node_id];
        std::pair<EdgeDescriptor, bool> new_edge = boost::add_edge(control_vertex, cntrl_c, controls_graph);
        controls_graph[new_edge.first].distance = distance;
        
    }
}


class dfs_interpolate_visitor :public boost::default_dfs_visitor
{

public:
	dfs_interpolate_visitor(Map_Double_SPtr _dem, Map_Double_SPtr _changes, boost::progress_display & _display_progress)
    :  down_cntrl_id(-1),
       up_cntrl_id(-1),
       dist_from_down_cntrl(-1),
       dist_from_up_cntrl(-1),
       dem(_dem),
       changes(_changes),
       display_progress(_display_progress)
	{
	}

	template < typename Vertex, typename Graph >
	void discover_vertex(Vertex u, Graph & g)
	{
        
        //Check if junction. If so, make a junction control node.
        typename boost::graph_traits<Graph>::degree_size_type num_out_edges = boost::out_degree(u, g);
        if (num_out_edges > 1)
        {
            g[u].type = JUNCT_CNTRL;
        }
        
        //If control, restart count and set new downstream control node.
        if (g[u].type > 1)
        {
            this->down_cntrl_id = u;
            this->dist_from_down_cntrl = 0;
        }
        //Else not junction, continue with distance count
        else
        {
            ++dist_from_down_cntrl;
        }
        
        //If end (no outlinks), then check is source node, reset up control id, and update up control distance accordingly.
        if (num_out_edges < 1)
        {
            this->up_cntrl_id = u; g[u].up_cntrl_id = u;
            this->dist_from_up_cntrl = 0; g[u].up_cntrl_dist = 0;
        }
        
        //Set node properties
        g[u].down_cntrl_dist = this->dist_from_down_cntrl;
        g[u].down_cntrl_id = this->down_cntrl_id;
        
	}
	
	template < typename Vertex, typename Graph >
	void finish_vertex(Vertex u, Graph & g)
	{
        //If control, restart count and set new downstream control node.
        if (g[u].type > 1)
        {
            this->up_cntrl_id = u;
            this->dist_from_up_cntrl = 0;
        }
        //Else not junction, continue with distance count
        else
        {
            ++dist_from_up_cntrl;
        }
        
        //Check if outflow node. If so, make a outflow control node.
        typename boost::graph_traits<Graph>::degree_size_type num_in_edges = boost::in_degree(u, g);
        if (num_in_edges < 1)
        {
            this->up_cntrl_id = u; g[u].up_cntrl_id = u;
            this->dist_from_up_cntrl = 0; g[u].up_cntrl_dist = 0;
        }


        
//		double elevation = dem->Get(g[u].row, g[u].col);
//		if (elevation > current_channel_elevation)
//		{
//			//then flow is wrong..., change to current_channel_elevation
//			//Calculate change
//			double change_val =  current_channel_elevation - elevation;
//			if (change_val < -300)
//			{
//				std::cout << "Something seems wrong here" << std::endl;
//			}
//			changes->Get(g[u].first, g[u].second) = change_val;
//			dem->Get(g[u].first, g[u].second) = current_channel_elevation;
//		}
//		current_channel_elevation = dem->Get(g[u].first, g[u].second);
//		if (current_channel_elevation < -300)
//		{
//			std::cout << "Something seems wrong here" << std::endl;
//		}
		++display_progress;
	}
	
private:
    int down_cntrl_id;
    int up_cntrl_id;
    int dist_from_down_cntrl;
    int dist_from_up_cntrl;
    
	Map_Double_SPtr dem;
	Map_Double_SPtr changes;
	boost::progress_display & display_progress;
};


#endif //CREEK_INCISOR_VISITOR
