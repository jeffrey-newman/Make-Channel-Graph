#include <boost/foreach.hpp>
#include "AddAdjacent.h"
#include "Neighbourhood.h"

void
SetControl(int row, int col, Map_Int_SPtr feature_map, ChannelNode & node_properties)
{
    int32_t val = feature_map->Get(row, col);
    if (val != feature_map->NoDataValue())
    {
        //channel_pixels.insert(std::make_pair(val, Position(i, j)));
        if (val == TERMINAL_CNTRL)
        {
            node_properties.type = TERMINAL_CNTRL;
        }
//        if (val == OUTFLW_CNTRL)
//        {
//            node_properties.type = OUTFLW_CNTRL;
//        }
    }

}

//void
//SetControl(int row, int col, Map_Int_SPtr feature_map, ChannelNode & node_properties)
//{
//    int32_t val = feature_map->Get(row, col);
//    if (val != feature_map->NoDataValue())
//    {
//        //channel_pixels.insert(std::make_pair(val, Position(i, j)));
//        if (val == OUTFLW_CNTRL)
//        {
//            node_properties.type = OUTFLW_CNTRL;
//        }
//        if (val == OUTFLW_CNTRL)
//        {
//            node_properties.type = OUTFLW_CNTRL;
//        }
//    }
//    
//}


void
SetIfSourceControl(VertexDescriptor v, Graph & g, Map_Int_SPtr feature_map)
{
    int32_t val = feature_map->Get(g[v].row, g[v].col);
    if (val == GUAGE_CNTRL)
    {
        g[v].type = GUAGE_CNTRL;
    }
}

void
SetIfSourceControl(VertexDescriptor v, Graph & g, Map_Int_SPtr feature_map, Map_Double_SPtr dem_map)
{
    //    // Search the neighbourhood.
    //    //Find out how many unprocessed adjacent creek features there are.
    //    ChannelNode & v_location = g[v];
    //    boost::shared_ptr<Set> Af = find_adjacents(feature_map, v_location.first, v_location.second, 1);
    //    if (Af->size() < 2)
    //    {
    //Then is source or sink, as only one adjacent cell.
    //Determine source or sink by gradient across nearest pixels
    int linkage_depth = 5;
    ChannelNode & previousNode = g[v];
    boost::shared_ptr<Set> tempAf = find_adjacents(feature_map, previousNode.row, previousNode.col, 1);
    ChannelNode & currentNode = tempAf->front();
    int source_count = 0;
    int sink_count = 0;
    
    while (linkage_depth > 0)
    {
        if (dem_map->Get(previousNode.row, previousNode.col) - dem_map->Get(currentNode.row, currentNode.col) > 0 )
        {
            // Then source
            ++source_count;
        }
        if (dem_map->Get(previousNode.row, previousNode.col) - dem_map->Get(currentNode.row, currentNode.col) > 0 )
        {
            // Then sink
            ++sink_count;
        }
        
        tempAf = find_adjacents(feature_map, currentNode.row, currentNode.col, 1);
        previousNode = currentNode;
        currentNode = tempAf->front();
        if (tempAf->size() > 2) linkage_depth = 0;
        else --linkage_depth;
    }
    
}


boost::optional<VertexDescriptor>
add_adjacent(VertexDescriptor v, Graph & g, Map_Double_SPtr dem_map, Map_Int_SPtr feature_map, Map_Bool_SPtr processed_map, std::list<ChannelNode> & outflow_pixels)
{
	// Search the neighbourhood.
	//Find out how many unprocessed adjacent creek features there are.
	ChannelNode & v_location = g[v];
	boost::shared_ptr<Set> Af = find_adjacents(feature_map, v_location.row, v_location.col, 1);
	std::vector<ChannelNode> adjacents;
	BOOST_FOREACH(ChannelNode & n_location, *Af)
	{
		// Test if adjacent pixel is processed and is creek
		double val = feature_map->Get(n_location.row, n_location.col);
		bool is_processed = processed_map->Get(n_location.row, n_location.col);
		//If not processed and a creek, add to graph.
		if ((val != feature_map->NoDataValue()) && !is_processed)
		{
			adjacents.push_back(n_location);
		}

	}
	//If there are no adjacents, then we are at the end of the channel
	if (adjacents.empty())
	{
//        v_location.type = OUTFLW_CNTRL;
        v_location.type = TERMINAL_CNTRL;
        outflow_pixels.push_back(g[v]);
		return (boost::none);
	}

	//If there are adjacents, add the highest adjacent as the next downstream node, and form a link.
	//Find highest adjacent, and add this as the downstream node
	ChannelNode highest;
	double highest_val = -999;
	BOOST_FOREACH(ChannelNode & location, adjacents)
	{
		double val = dem_map->Get(location.row, location.col);
		if (val > highest_val)
		{
			highest = location;
			highest_val = val;
		}
	}
	//Add vertex for node
	VertexDescriptor new_v = boost::add_vertex(g);
    g[new_v] = highest;
    
    //Check if guage control point and add this info to the graph.
    int32_t val = feature_map->Get(highest.row, highest.col);
    //channel_pixels.insert(std::make_pair(val, Position(i, j)));
    if (val == GUAGE_CNTRL)
    {
        g[new_v].type = GUAGE_CNTRL;
    }
    
    
    
    
	processed_map->Get(highest.row, highest.col) = true; // We have now processed this cell.
	// Add link joining to source node.
	EdgeDescriptor edgeD;
	bool inserted;
    //We want the inverse direction.
	//boost::tie(edgeD, inserted) = boost::add_edge(v, new_v, g);
    boost::tie(edgeD, inserted) = boost::add_edge(new_v, v, g);
	if (inserted)
	{
		return (new_v);
	}
	else
	{
		return (boost::none);
	}
}
