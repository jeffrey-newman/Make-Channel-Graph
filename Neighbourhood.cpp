//
//  neighbourhood.cpp
//  move_creek
//
//  Created by a1091793 on 23/01/2015.
//  Copyright (c) 2015 University of Adelaide. All rights reserved.
//

#include "Neighbourhood.h"
#include <stdio.h>
#include <cmath>
#include <algorithm>



template <typename DataType>
boost::shared_ptr<Set> get_neighbourhood(boost::shared_ptr<Map_Matrix<DataType> > map, unsigned long i, unsigned long j, unsigned long size)
{
    int ii = (int) i;
    int jj = (int) j;
    int sizei = (int) size;
    boost::shared_ptr<Set> neighbourhood(new Set);
    DataType noData = map->NoDataValue();
    
    int end_y = (int) std::min(ii+sizei, (int) (map->NRows() - 1) );
    int begin_y = (int) std::max(ii-sizei, 0);
    
    int end_x = (int) std::min(jj+sizei, (int) (map->NCols() - 1) );
    int begin_x = (int) std::max(jj-sizei, 0);
    
    for (int ti = begin_y; ti <= end_y; ++ti)
    {
        for (int tj = begin_x; tj <= end_x; ++tj)
        {
            if (euclidean_distance(ti, tj, ii, jj) < double(sizei) && map->Get(ti, tj) != noData )
            {
                ChannelNode node;
                node.row = ti;
                node.col = tj;
                neighbourhood->push_back(node);
            }
        }
    }
    return (neighbourhood);
}

template <typename DataType>
boost::shared_ptr<Set> find_adjacents(boost::shared_ptr<Map_Matrix<DataType> > map, unsigned long i, unsigned long j, unsigned long size)
{
    int ii = (int) i;
    int jj = (int) j;
    int radius = (int) size;
    boost::shared_ptr<Set> adjacents(new Set);
    DataType noData = map->NoDataValue();
    
	int end_y = (int)std::min(ii + radius, (int)(map->NRows() - 1));
	int begin_y = (int)std::max(ii - radius, 0);
    
	int end_x = (int)std::min(jj + radius, (int)(map->NCols() - 1));
	int begin_x = (int)std::max(jj - radius, 0);
    
    for (int ti = begin_y; ti <= end_y; ++ti)
    {
        for (int tj = begin_x; tj <= end_x; ++tj)
        {
			DataType val = (DataType) map->Get(ti, tj);
			//double dist = euclidean_distance(ti, tj, ii, jj);
			//if (dist < (double) radius && val != noData) adjacents->push_back(std::make_pair(ti - ii, tj - jj));
            ChannelNode node;
            node.row = ti;
            node.col = tj;
			if (val != noData) adjacents->push_back(node);
        }
    }
    return (adjacents);
}
