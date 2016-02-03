#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <map>
#include <list>
#include <limits>  
#include <cmath>

#include <GDAL/gdal.h>

#include <boost/filesystem.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/graph/undirected_dfs.hpp>


#include <boost/foreach.hpp>
#include <boost/progress.hpp>


#include "Types.h"
#include "ReadInMap.h"
#include "Print.h"
#include "AddAdjacent.h"
#include "Neighbourhood.h"
#include "GraphSearches.h"
#include "PrintGraphsToFile.h"



int main(int argc, char **argv)
{


	/**********************************/
	/*        Program options         */
	/**********************************/
	// Need to specify elevation grid
	// Need to specify channel 
	std::string feature_file;
	std::string dem_file;
	double gap_threshold;
	std::string out_file_name;
	std::string control_file_name;

    unsigned int trim_level;

	namespace prog_opt = boost::program_options;
	prog_opt::options_description desc("Allowed options");
	desc.add_options()
		("help,h", "produce help message")
		("feature-map,f", prog_opt::value<std::string>(&feature_file), "path of the gdal capatible raster feature file")
		("dem-map,d", prog_opt::value<std::string>(&dem_file), "path of the gdal capatible elevation data file")
		("trim-branches,t", prog_opt::value<unsigned int>(&trim_level)->default_value(0), "remove branches if they are composed with a number of pixels less than this amount")
		("bridge-gaps-threshold,b", prog_opt::value<double>(&gap_threshold)->default_value(0), "Threshold distance in which channel graph will bridge a gap in the channel raster")
		("graph-file,g", prog_opt::value<std::string>(&out_file_name)->default_value("channel_graph"), "Name of channel graph file that will be produced")
		("terminal-nodes,n", prog_opt::value<std::string>(&control_file_name)->default_value("terminal_nodes.txt"), "Name of file in which the terminal nodes present in graph will be written to");

	prog_opt::variables_map vm;
	prog_opt::store(prog_opt::parse_command_line(argc, argv, desc), vm);
	prog_opt::notify(vm);
	if (vm.count("help")) 
	{
		std::cout << desc << "\n";
		return 1;
	}


	boost::filesystem::path feature_file_path(feature_file);
	boost::filesystem::path dem_file_path(dem_file);
	


	// Check file exists
	if (!fs::exists(feature_file_path))
	{
		std::stringstream ss;
		ss << feature_file_path << " does not exist";
		throw std::runtime_error(ss.str());
		return (EXIT_FAILURE);
	}

	if (!fs::exists(dem_file_path))
	{
		std::stringstream ss;
		ss << dem_file_path << " does not exist";
		throw std::runtime_error(ss.str());
		return (EXIT_FAILURE);
	}
 

    /**********************************/
    /*       Create graph object      */
    /**********************************/
    Graph channel_grph;
    
    

	/**********************************/
	/*         Read in maps           */
	/**********************************/
	std::cout << "\n\n*************************************\n";
	std::cout <<     "*             Read in maps          *\n";
	std::cout <<     "*************************************" << std::endl;
	//First the channel
	std::tuple<boost::shared_ptr<Map_Matrix<int32_t> >, std::string, GeoTransform> gdal_feature_map = read_in_map<int32_t>(feature_file_path, GDT_Int32, NO_CATEGORISATION);
	boost::shared_ptr<Map_Matrix<int32_t> > feature_map(std::get<0>(gdal_feature_map));
	std::string & featureWKTprojection(std::get<1>(gdal_feature_map));
	GeoTransform & featureTransform(std::get<2>(gdal_feature_map));

	// Second the elevation data
	std::tuple<boost::shared_ptr<Map_Matrix<double> >, std::string, GeoTransform> gdal_dem_map = read_in_map<double>(dem_file_path, GDT_Float64, NO_CATEGORISATION);
	boost::shared_ptr<Map_Matrix<double> > dem_map(std::get<0>(gdal_dem_map));
	std::string & demWKTprojection(std::get<1>(gdal_dem_map));
	GeoTransform & demTransform(std::get<2>(gdal_dem_map));

	//Check maps for consistency (same dimensions)
	if (feature_map->NCols() != dem_map->NCols())
	{
		throw std::runtime_error("Number of columns in the two comparison maps non-identical");
	}

	if (feature_map->NRows() != dem_map->NRows())
	{
		throw std::runtime_error("Number of rows in the two comparison maps non-identical");
	}

	/**********************************/
	/*  Create list of channel pixels */
	/**********************************/
	std::cout << "\n\n*************************************\n";
	std::cout <<     "*  Making list of channel pixels:   *\n";
	std::cout <<     "*************************************" << std::endl;
	boost::progress_display show_progress1((feature_map->NRows()*feature_map->NCols()));
	//std::multimap<double, Position> channel_pixels;
    
    std::map<int, std::map<int, VertexDescriptor>  > channel_pixels;
//    bool hasNoData = feature_map->HasNoDataValue();
    int32_t noDataVal = feature_map->NoDataValue();
    
    int node_id_count = 0;
    for (unsigned int i = 0; i < feature_map->NRows(); ++i)
    {
        for (unsigned int j = 0; j < feature_map->NCols(); ++j)
        {
            int32_t val = feature_map->Get(i, j);
            if (val != noDataVal)
            {
                //channel_pixels.insert(std::make_pair(val, Position(i, j)));
                VertexDescriptor v = boost::add_vertex(channel_grph);
                channel_grph[v].row = i;
                channel_grph[v].col = j;
                channel_grph[v].node_id = node_id_count++;
                channel_grph[v].elevation = dem_map->Get(i,j);
//                if (val == GUAGE_CNTRL)
//                {
//                    channel_grph[v].type = GUAGE_CNTRL;
//                }
                channel_pixels[i].insert(std::make_pair(j, v));
            }
            ++show_progress1;
        }
    }



	//// At this point we no longer need the data in the dem_map and can safely delete this data from memory
	dem_map.reset(new Map_Double);
	std::get<0>(gdal_dem_map).reset(new Map_Double);
	

	/**********************************/
	/*    Create graph of channel     */
	/**********************************/
	std::cout << "\n\n*************************************\n";
	std::cout <<     "*       Making channel graph:       *\n";
	std::cout <<     "*************************************" << std::endl;
	
    

    // Add edges between all adjacent channel cells
	boost::progress_display show_progress2(channel_pixels.size());
    
    typedef std::pair<const int, std::map<int, VertexDescriptor> > RowIteratorType;
    typedef std::pair<const int, VertexDescriptor> ColIteratorType;
    BOOST_FOREACH(RowIteratorType & row_it, channel_pixels)
    {
        BOOST_FOREACH(ColIteratorType & channel_vertex, row_it.second)
        {
            ChannelNode & v_location = channel_grph[channel_vertex.second];
            boost::shared_ptr<Set> Af = find_adjacents(feature_map, v_location.row, v_location.col, 1);
            BOOST_FOREACH(ChannelNode & n_location, *Af)
            {
                int row = n_location.row;
                int col = n_location.col;
                VertexDescriptor adjacent_vertex = channel_pixels[row][col];
                boost::add_edge(channel_vertex.second, adjacent_vertex, channel_grph);
            }
        }
        ++show_progress2;
    }

	//// At this point we no longer need the feature map data and can safely delete it from memory
	std::get<0>(gdal_feature_map).reset(new Map_Int);
	feature_map.reset(new Map_Int);


    
    /**********************************/
    /*     Remove loops in graph      */
    /**********************************/
	std::cout << "\n\n*************************************\n";
	std::cout <<     "*      Removing loops in graph      *\n";
	std::cout <<     "*************************************" << std::endl;
    VertexDescMap idxMap;
    VertexIDMap idMap;
    boost::associative_property_map<VertexDescMap> indexMap(idxMap);
    VertexIterator di, dj;
    boost::tie(di, dj) = boost::vertices(channel_grph);
    for(int i = 0; di != dj; ++di,++i)
    {
        boost::put(indexMap, (*di), channel_grph[*di].node_id);
        idMap.insert(std::make_pair(channel_grph[*di].node_id, (*di)));
    }
    
    typedef std::map<VertexDescriptor, boost::default_color_type> VertexColourMap;
    typedef std::map<EdgeDescriptor, boost::default_color_type> EdgeColourMap;
    VertexColourMap vertex_colour_map;
    EdgeColourMap edge_colour_map;
    boost::associative_property_map<EdgeColourMap> edge_colours(edge_colour_map);
    boost::associative_property_map<VertexColourMap> vertex_colours(vertex_colour_map);
    

    std::vector<EdgeDescriptor> removal_list;
    CycleTerminator<EdgeDescriptor, Graph> vis(removal_list);
    boost::undirected_dfs(channel_grph,boost::visitor(vis).vertex_color_map(vertex_colours).edge_color_map(edge_colours).root_vertex(*vertices(channel_grph).first).vertex_index_map(indexMap));
    BOOST_FOREACH(EdgeDescriptor edge, removal_list)
    {
        boost::remove_edge(edge, channel_grph);
    }

    /**********************************/
    /*   Identify controls in Graph   */
    /**********************************/
    std::cout << "\n\n*************************************\n";
    std::cout <<     "*  Identifying controls in graph:   *\n";
    std::cout <<     "*************************************" << std::endl;
    
    //Identify terminal nodels and junctional nodes based on number of out_edges.
	std::list<VertexDescriptor> terminal_nodes;
    std::pair<VertexIterator, VertexIterator> vp;
    for (vp = boost::vertices(channel_grph); vp.first != vp.second; ++vp.first)
    {
        OutDegreeType num_edges = boost::out_degree(*(vp.first), channel_grph);
        if (num_edges == 1)
        {
            channel_grph[*vp.first].type = TERMINAL_CNTRL;
			terminal_nodes.push_back(*vp.first);
        }
        if (num_edges > 2)
        {
            channel_grph[*vp.first].type = JUNCT_CNTRL;
        }
    }

	/**********************************/
	/*   Bridging gaps in graph:   */
	/**********************************/
	std::cout << "\n\n*************************************\n";
	std::cout << "*  Bridging gaps in graph:   *\n";
	std::cout << "*************************************" << std::endl;
	if (gap_threshold > 0)
	{
		for (std::list<VertexDescriptor>::const_iterator it = terminal_nodes.begin(); it != terminal_nodes.end(); ++it)
		{
			double min_dist = std::numeric_limits<double>::max();
			VertexDescriptor adjacent_v = *it;

			for (std::list<VertexDescriptor>::const_iterator it2 = terminal_nodes.begin(); it2 != terminal_nodes.end(); ++it2)
			{
				//distance between it and it2
				if (channel_grph[*it].node_id != channel_grph[*it].node_id)
				{
					int row1 = channel_grph[*it].row;
					int col1 = channel_grph[*it].col;
					int row2 = channel_grph[*it2].row;
					int col2 = channel_grph[*it2].col;
					double dist = sqrt(pow(double(row2 - row1), 2) + pow(double(col2 - col1), 2));
					if (dist < min_dist)
					{
						min_dist = dist;
						adjacent_v = *it2;
					}
				}
			}

			if (min_dist <= gap_threshold)
			{
				boost::add_edge(*it, adjacent_v, channel_grph);
				channel_grph[*it].type = NOT_CNTRL;
				channel_grph[adjacent_v].type = NOT_CNTRL;
			}
		}
	}
    
    /*****************************************************/
    /*   Remove leaf tributaries if length < trim_level  */
    /*****************************************************/
	std::cout << "\n\n******************************************************\n";
	std::cout <<     "*   Remove leaf tributaries if length < trim_level   *\n";
	std::cout <<     "******************************************************" << std::endl;
    for (int it_trim_levl = 0; it_trim_levl <= trim_level ; ++it_trim_levl)
    {
		std::set<EdgeDescriptor> edge_removal_list;
		std::set<VertexDescriptor> vertex_removal_list;
        for (vp = boost::vertices(channel_grph); vp.first != vp.second; ++vp.first)
        {
            if(channel_grph[*vp.first].type == TERMINAL_CNTRL)
            {
                leafBranchSearch(channel_grph, *vp.first, it_trim_levl, edge_removal_list, vertex_removal_list);
            }
        }
        BOOST_FOREACH(EdgeDescriptor e, edge_removal_list)
        {
            //std::cout << "deleting edge ( " << channel_grph[boost::source(e,channel_grph)].node_id << " , " << channel_grph[boost::target(e,channel_grph)].node_id << " )" << std::endl;
            boost::remove_edge(e, channel_grph);
        }
        BOOST_FOREACH(VertexDescriptor v, vertex_removal_list)
        {
            //std::cout << "deleting vertex " << v << std::endl;
            boost::remove_vertex(v, channel_grph);
        }
        for (vp = boost::vertices(channel_grph); vp.first != vp.second; ++vp.first)
        {
            if(channel_grph[*vp.first].type == JUNCT_CNTRL)
            {
                if (boost::out_degree(*vp.first, channel_grph) < 3)
                {
                    channel_grph[*vp.first].type = NOT_CNTRL;
                }
            }
            if (boost::out_degree(*(vp.first), channel_grph) == 1)
            {
                channel_grph[*vp.first].type = TERMINAL_CNTRL;
            }
        }
    }

    


//    /********************************************/
//    /*        Create control node graph         */
//    /********************************************/
//    std::cout << "\n\n*************************************\n";
//    std::cout <<     "*    Creating control node graph    *\n";
//    std::cout <<     "*************************************" << std::endl;
//    
//    Graph controls_grph;
//    VertexIDMap node_id_map;
////    std::pair<VertexIterator, VertexIterator> vp;
//    for (vp = boost::vertices(channel_grph); vp.first != vp.second; ++vp.first)
//    {
//        if (channel_grph[*(vp.first)].type > 1)
//        {
//            VertexDescriptor v = boost::add_vertex(controls_grph);
//            controls_grph[v] = channel_grph[*(vp.first)];
//            node_id_map.insert(std::make_pair(channel_grph[*(vp.first)].node_id, v));
//            channel_grph[*(vp.first)].up_cntrl_id = channel_grph[*(vp.first)].node_id;
//            channel_grph[*(vp.first)].up_cntrl_dist = 1;
//            channel_grph[*(vp.first)].down_cntrl_id = channel_grph[*(vp.first)].node_id;
//            channel_grph[*(vp.first)].down_cntrl_dist = 1;
//        }
//    }
//    
//    // Create links in controls graph using a search across the channels graph at each point.
//    for (vp = boost::vertices(controls_grph); vp.first != vp.second; ++vp.first)
//    {
//        int control_id = controls_grph[*(vp.first)].node_id;
//        VertexDescriptor channel_vertex = idMap[control_id];
//        controlSearch(channel_grph, channel_vertex, controls_grph, *(vp.first), node_id_map);
//    }
    

    
    /********************************************/
    /*       Print Control graphs to file       */
    /********************************************/
    //std::cout << "\n\n*************************************\n";
    //std::cout <<     "*  Printing control Graph to file   *\n";
    //std::cout <<     "*************************************" << std::endl;
    //std::string controls_file_name = "control_graph";
    //printGraphsToFile(controls_grph, controls_file_name);
    
    /********************************************/
    /*       Print Channel graphs to file       */
    /********************************************/
    std::cout << "\n\n*************************************\n";
    std::cout <<     "*  Printing channel Graph to file   *\n";
    std::cout <<     "*************************************" << std::endl;
	printGraphsToFile(channel_grph, out_file_name);
    
	/********************************************/
	/*    Print terminal node ids to file       */
	/********************************************/
	std::cout << "\n\n*************************************\n";
	std::cout <<     "*   Print terminal node ids to file *\n";
	std::cout <<     "*************************************" << std::endl;
	std::ofstream tn_fs(control_file_name.c_str());
	if (tn_fs.is_open())
	{
		for (vp = boost::vertices(channel_grph); vp.first != vp.second; ++vp.first)
		{
			if (channel_grph[*(vp.first)].type == TERMINAL_CNTRL)
			{
				tn_fs << channel_grph[*(vp.first)].node_id << "\tterminal_type\t3\t" << channel_grph[*(vp.first)].elevation << "\n";
			}
		}
	}

    return (EXIT_SUCCESS);
}

	