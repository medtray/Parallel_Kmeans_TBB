#include <iostream>
#include <vector>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <algorithm>
#include <chrono>
#include <fstream>
#include <sstream>
#include "tbb/tbb.h"
#include "tbb/task_scheduler_init.h"
#include "tbb/enumerable_thread_specific.h"

using namespace std;
using namespace tbb;





class Point
{
private:
	int id_point, id_cluster;
	vector<double> values;
	int total_values;
	string name;

public:
	Point(int id_point, vector<double>& values, string name = "")
	{
		this->id_point = id_point;
		total_values = values.size();

		for(int i = 0; i < total_values; i++)
			this->values.push_back(values[i]);

		this->name = name;
		id_cluster = -1;
	}

	int getID()
	{
		return id_point;
	}

	vector<double> getValues()
	{
		return values;
	}

	void setCluster(int id_cluster)
	{
		this->id_cluster = id_cluster;
	}

	int getCluster()
	{
		return id_cluster;
	}

	double getValue(int index)
	{
		return values[index];
	}

    void updateValue(int index,double val)
    {
        values[index]=val;
    }

	void updateValues(vector<double> val)
    {
        values.assign(val.begin(),val.end()); 
    }

	int getTotalValues()
	{
		return total_values;
	}

	void addValue(double value)
	{
		values.push_back(value);
	}

	string getName()
	{
		return name;
	}
};


class sum_and_count
{
    public:
    size_t count;
    Point sum;
    Point init_point;
    sum_and_count(Point p):count(0),sum(p),init_point(p)
        
    {
    }

    public:

    void tally(Point p ) {
         for (int j=0;j<sum.getTotalValues();j++)
        {
            sum.updateValue(j,sum.getValue(j)+p.getValue(j));
        } 

		/* vector<double> inter=sum.getValues();
		
		 std::transform (inter.begin(), inter.end(), p.getValues().begin(), inter.begin(), std::plus<double>());
		 sum.updateValues(inter); */
        ++count;
        }

        void clear() {
            sum = init_point;
            count = 0;
        }

        Point mean(){
            for (int j=0;j<sum.getTotalValues();j++)
            {
                sum.updateValue(j,sum.getValue(j)/count);
            }
            return sum;
        }


    void operator+=(sum_and_count &other ) {
        for (int j=0;j<sum.getTotalValues();j++)
        {
            sum.updateValue(j,sum.getValue(j)+other.sum.getValue(j));
        }
        count += other.count;
    }

};


class view {
        view( const view& v );
        void operator=( const view& v ); 
        public:
        std::vector<sum_and_count> array;
        size_t change;
        view(std::vector<sum_and_count> array2) :array(array2), change(0) {}
        ~view() {array.clear();}
        
        };

typedef tbb::enumerable_thread_specific<view> tls_type;

void reduce_local_counts_to_global_count( tls_type& tls, view& global ) {
    global.change = 0;
    for( auto i=tls.begin(); i!=tls.end(); ++i ) {
        view& v = *i;
        global.change += v.change;
        v.change = 0;
    }
}

void reduce_local_sums_to_global_sum( size_t k, tls_type& tls, view& global ) {
for( auto i=tls.begin(); i!=tls.end(); ++i ) {
    view& v = *i;
    for( size_t j=0; j<k; ++j ) {
        global.array[j] += v.array[j];
        
        v.array[j].clear();
    }
    }
}

class Cluster
{
private:
	int id_cluster;
	vector<double> central_values;
	vector<Point> points;

public:
	Cluster(int id_cluster, Point point)
	{
		this->id_cluster = id_cluster;

		int total_values = point.getTotalValues();

		for(int i = 0; i < total_values; i++)
			central_values.push_back(point.getValue(i));

		points.push_back(point);
	}

	void addPoint(Point point)
	{
		points.push_back(point);
	}

	bool removePoint(int id_point)
	{
		int total_points = points.size();

		for(int i = 0; i < total_points; i++)
		{
			if(points[i].getID() == id_point)
			{
				points.erase(points.begin() + i);
				return true;
			}
		}
		return false;
	}

	double getCentralValue(int index)
	{
		return central_values[index];
	}

	void setCentralValue(int index, double value)
	{
		central_values[index] = value;
	}

	Point getPoint(int index)
	{
		return points[index];
	}

	int getTotalPoints()
	{
		return points.size();
	}

	int getID()
	{
		return id_cluster;
	}
};

class KMeans
{
private:
	int K; // number of clusters
	int total_values, total_points, max_iterations,n_threads,grain_size;
	vector<Cluster> clusters;

	// return ID of nearest center (uses euclidean distance)
	int getIDNearestCenter(Point point)
	{
		double sum = 0.0, min_dist;
		int id_cluster_center = 0;

		for(int i = 0; i < total_values; i++)
		{
			sum += pow(clusters[0].getCentralValue(i) -
					   point.getValue(i), 2.0);
		}

		min_dist = sqrt(sum);

		for(int i = 1; i < K; i++)
		{
			double dist;
			sum = 0.0;

			for(int j = 0; j < total_values; j++)
			{
				sum += pow(clusters[i].getCentralValue(j) -
						   point.getValue(j), 2.0);
			}

			dist = sqrt(sum);

			if(dist < min_dist)
			{
				min_dist = dist;
				id_cluster_center = i;
			}
		}

		return id_cluster_center;
	}

public:
	KMeans(int K, int total_points, int total_values, int max_iterations,int n_threads, int grain_size)
	{
		this->K = K;
		this->total_points = total_points;
		this->total_values = total_values;
		this->max_iterations = max_iterations;
		this->n_threads = n_threads;
		this->grain_size = grain_size;
	}

	void run(vector<Point> & points)
	{

        vector<double> init_values;

        for(int j = 0; j < total_values; j++)
        {
            init_values.push_back(0);
        }

        Point init_point(0, init_values);

        std::vector<sum_and_count> ss;


        for (int i=0;i<K;i++)
            ss.push_back(sum_and_count(init_point));



        
        tls_type tls([&]{return ss;});
        view global(ss);

        
        auto begin = chrono::high_resolution_clock::now();
        
		if(K > total_points)
			return;

		vector<int> prohibited_indexes;
		tbb::task_scheduler_init init(n_threads);
		// choose K distinct values for the centers of the clusters
		int index_point;
		for(int i = 0; i < K; i++)
		{
			while(true)
			{
				//int index_point = rand() % total_points;
				index_point = K*i+1;

				if(find(prohibited_indexes.begin(), prohibited_indexes.end(),
						index_point) == prohibited_indexes.end())
				{
					prohibited_indexes.push_back(index_point);
					points[index_point].setCluster(i);
					Cluster cluster(i, points[index_point]);
					clusters.push_back(cluster);
					break;
				}
			}
		}
        auto end_phase1 = chrono::high_resolution_clock::now();
        
		int iter = 1;
		if (grain_size==0)
			grain_size=total_points/n_threads;

		cout << "grain_size= " << grain_size << "\n\n";

		while(true)
		{
			bool done = true;


                        // Compute new clusters and their local sums
            tbb::parallel_for(
            tbb::blocked_range<size_t>(0,total_points,grain_size),
            [&]( tbb::blocked_range<size_t> r ) {
                view& v = tls.local();
                for( size_t i=r.begin(); i!=r.end(); ++i ) {
                    // “Reassign step”: Find index of centroid closest to points [i]
                    int id_nearest_center = getIDNearestCenter(points[i]);
                    int id_old_cluster = points[i].getCluster();
                    if( id_nearest_center!=id_old_cluster) {
                    points[i].setCluster(id_nearest_center);
                    ++v.change;
                    }
                    // “Sum step”
                    
                    v.array[id_nearest_center].tally(points[i]);
                }
            }

            ,tbb::simple_partitioner());

            // Reduce local counts to global count
            reduce_local_counts_to_global_count(tls, global);

            if(global.change!=0)
            {
                // Reduce local sums to global sum
                reduce_local_sums_to_global_sum(K, tls, global);
                // “Divide step”: Compute centroids from global sums

                for(int i = 0; i < K; i++)
			{
                Point center=global.array[i].mean();
				for(int j = 0; j < total_values; j++)
				{
                    clusters[i].setCentralValue(j, center.getValue(j));
                }
                global.array[i].clear();
            }

            
                done = false;
            }

			if(done == true || iter >= max_iterations)
			{
				cout << "Break in iteration " << iter << "\n\n";
				break;
			}

			iter++;
		}
        auto end = chrono::high_resolution_clock::now();

		// shows elements of clusters
		for(int i = 0; i < K; i++)
		{
			int total_points_cluster =  clusters[i].getTotalPoints();

			cout << "Cluster " << clusters[i].getID() + 1 << endl;
			for(int j = 0; j < total_points_cluster; j++)
			{
				cout << "Point " << clusters[i].getPoint(j).getID() + 1 << ": ";
				for(int p = 0; p < total_values; p++)
					cout << clusters[i].getPoint(j).getValue(p) << " ";

				string point_name = clusters[i].getPoint(j).getName();

				if(point_name != "")
					cout << "- " << point_name;

				cout << endl;
			}

			cout << "Cluster values: ";

			for(int j = 0; j < total_values; j++)
				cout << clusters[i].getCentralValue(j) << " ";

			cout << "\n\n";
            cout << "TOTAL EXECUTION TIME = "<<std::chrono::duration_cast<std::chrono::microseconds>(end-begin).count()<<"\n";
            
            cout << "TIME PHASE 1 = "<<std::chrono::duration_cast<std::chrono::microseconds>(end_phase1-begin).count()<<"\n";
            
            cout << "TIME PHASE 2 = "<<std::chrono::duration_cast<std::chrono::microseconds>(end-end_phase1).count()<<"\n";
		}
	}
};

/** Print some helpful usage information */
void usage() {
  using std::cout;
  cout << "Parallel K-Means\n";
  cout << "  Usage: gauss [options]\n";
  cout << "    -t       : number of threads (default 4)\n";
  cout << "    -g       : grain size (default total_points/n_threads)\n";
  cout << "    -h       : print this message\n";
}

int main(int argc, char *argv[])
{

	int NThreads=4;
  int grain_size=0;

  // Parse the command line options:
  int o;
  while ((o = getopt(argc, argv, "t:g:h")) != -1) {
    switch (o) {
    case 't':
      NThreads = atoi(optarg);
      break;
    case 'g':
      grain_size = atoi(optarg);
      break;
    case 'h':
      usage();
	  exit(-1);
      break;
    default:
      usage();
      exit(-1);
    }
  }
	srand (time(NULL));

	int total_points, total_values, K, max_iterations, has_name;
	int line_in_file=-1;
	int index_point=-1;

	vector<Point> points;
	string point_name;

	string FILENAME="./datasets/dataset5.txt";

	std::ifstream file(FILENAME);
	if (file.is_open()) {
		std::string line;
		while (getline(file, line)) {
			line_in_file+=1;
			if (line_in_file==0)
			{
				vector <string> tokens; 
				std::stringstream check1(line); 
				string intermediate; 
				while(getline(check1, intermediate, ' ')) 
				{ 
					tokens.push_back(intermediate); 
				} 

				total_points=stoi(tokens[0]);
				total_values=stoi(tokens[1]);
				K=stoi(tokens[2]);
				max_iterations=stoi(tokens[3]);
				has_name=stoi(tokens[4]);

			}

			else
			{
				index_point+=1;
				vector <string> tokens; 
				std::stringstream check1(line); 
				string intermediate; 
				while(getline(check1, intermediate, ' ')) 
				{ 
					tokens.push_back(intermediate); 
				}

				vector<double> values;

				for(int j = 0; j < total_values; j++)
				{
					double value = atof(tokens[j].c_str());
					values.push_back(value);
				}

				if(has_name)
				{
					Point p(index_point, values, tokens[total_values]);
					points.push_back(p);
				}
				else
				{
					Point p(index_point, values);
					points.push_back(p);
				} 
				
			}
			
			
		}
		file.close();
	}

	KMeans kmeans(K, total_points, total_values, max_iterations,NThreads,grain_size);
	kmeans.run(points);

	return 0;
}
