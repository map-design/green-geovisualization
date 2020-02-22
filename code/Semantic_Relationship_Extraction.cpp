#include "stdafx.h"
#include <algorithm>

//input map color and cluster
void iniCluster(string path)
{
	ifstream infile;
	infile.open(path);

	string lineStr;
	vector<vector<string>> strArray;
	int colorid = 0;
	while (getline(infile, lineStr))
	{
		stringstream ss(lineStr);
		string str;
		vector<string> lineArray;

		while (getline(ss, str, ' '))
			lineArray.push_back(str);
		
		metacolor c;
		c.R = stoi(lineArray[0]);
		c.G= stoi(lineArray[1]);
		c.B = stoi(lineArray[2]);
		c.weight = stof(lineArray[3]);
		c.colorID = colorid++;
		calmetaColor(c);
		inputColor.emplace_back(c);
	}

	
	D = inputColor.size();
	calCluster(inputColor);
	

	for (int i = 0; i < colorCluster.size(); i++)
	{
		for (int j = 0; j < colorCluster[i].colorIDs.size(); j++)
		{
			inputColor[colorCluster[i].colorIDs[j]].clusterID = i;
		}
	}
	for (int i = 0; i < colorCluster.size(); i++)
	{
		if (colorCluster[i].type == 1)
		{
			vector<metacolor> orcolors;
			for (int j = 0; j < colorCluster[i].colorIDs.size(); j++)
			{
				orcolors.emplace_back(inputColor[colorCluster[i].colorIDs[j]]);
			}
			std::sort(orcolors.begin(), orcolors.end(), [](const metacolor& x, const metacolor& y)
			{
				return x.L > y.L;
			});
			
			inputColor[orcolors[0].colorID].selected = true;
			inputColor[orcolors[orcolors.size()-1].colorID].selected = true;
			inputColor[orcolors[orcolors.size() / 2].colorID].selected = true;
		}
	}

}

//Calculate color clusters
void calCluster(vector<metacolor>& colors)
{
	
	vector<point> m_lab;
	for (int i = 0; i < D; i++)
	{
		point temp;
		temp.id = i;
		temp.z = colors[i].L;
		temp.x = colors[i].A;
		temp.y = colors[i].Bright;
		temp.proportion = colors[i].weight;
		m_lab.push_back(temp);
	}

	DBSCAN(m_lab, DBSCANEps, 1, 1);

	map<int, int> clusters;

	for (int i = 0; i < D; i++)
	{
		cluster cc;
		m_lab[1].cluster = 0;
		m_lab[2].cluster = 0;

		if (m_lab[i].cluster == 0)
		{
			cc.type = 0;
			cc.colorIDs.emplace_back(i);
			colorCluster.emplace_back(cc);
		}
		else
		{
			map<int, int>::iterator it = clusters.find(m_lab[i].cluster);
			if (it != clusters.end())
			{
				colorCluster[(*it).second].colorIDs.emplace_back(i);
			}
			else
			{
				clusters[m_lab[i].cluster] = colorCluster.size();
				cc.colorIDs.emplace_back(i);
				colorCluster.emplace_back(cc);
			}
		}
	}
	for (int i = 0; i < colorCluster.size(); i++)
	{
		if (colorCluster[i].colorIDs.size() > 2)
			colorCluster[i].type = 1;
		if (colorCluster[i].colorIDs.size() == 2)
			colorCluster[i].type = 2;
	}
}

float squareDeviation(point a, point b)
{
	float m_c;
	m_c = (sqrt(pow((a.z - b.z), 2) + pow((a.x - b.x), 2) + pow((a.y - b.y), 2)));
	return m_c;
}

void DBSCAN(vector<point>& dataset, float Eps, int MinPts, int clustertime)
{
	int len = dataset.size();
	for (int i = 0; i < len; i++)
	{
		dataset[i].cluster = 0;
		dataset[i].pointType = 1;
		dataset[i].pts = 0;
		dataset[i].j = 0;
		dataset[i].visited = 0;
		dataset[i].corepts.clear();
	}
	//calculate pts
	for (int i = 0; i < len; i++)
	{
		for (int j = i + 1; j < len; j++)
		{
			if (squareDeviation(dataset[i], dataset[j]) < Eps)
			{
				dataset[i].pts++;
				dataset[j].pts++;
			}
		}
	}
	//core point 
	for (int i = 0; i < len; i++)
	{
		if (dataset[i].pts >= MinPts)
		{
			dataset[i].pointType = 3;
			corePoint.push_back(dataset[i]);
		}
	}
	//joint core point
	for (int i = 0; i < corePoint.size(); i++)
	{
		for (int j = i + 1; j < corePoint.size(); j++)
		{
			if (squareDeviation(corePoint[i], corePoint[j]) < Eps)
			{
				corePoint[i].corepts.push_back(j);
				corePoint[j].corepts.push_back(i);
			}
		}
	}
	int k = 1;
	for (int i = 0; i < corePoint.size(); i++)
	{
		stack<point*> ps;
		if (corePoint[i].visited == 1) continue;
		corePoint[i].cluster = k;
		ps.push(&corePoint[i]);
		point *v;
		while (!ps.empty())
		{
			v = ps.top();
			v->visited = 1;
			ps.pop();
			for (int j = 0; j < v->corepts.size(); j++)
			{
				if (corePoint[v->corepts[j]].visited == 1) continue;
				corePoint[v->corepts[j]].cluster = corePoint[i].cluster;
				corePoint[v->corepts[j]].visited = 1;
				ps.push(&corePoint[v->corepts[j]]);
			}
		}
		k = k + 1;
	}
	//border point,joint border point to core point
	for (int i = 0; i < len; i++)
	{
		if (dataset[i].pointType == 3) continue;
		for (int j = 0; j < corePoint.size(); j++)
		{
			if (squareDeviation(dataset[i], corePoint[j]) < Eps)
			{
				dataset[i].pointType = 2;
				dataset[i].cluster = corePoint[j].cluster;
				break;
			}
		}
	}
	//output
	fstream Fclustering;
	if (clustertime == 1)
	{
		Fclustering.open("result\\Fclustering.txt", ios::out);
		for (int i = 0; i < len; i++)
		{
			if (dataset[i].pointType == 2)
			{
				Fclustering << dataset[i].id << " " << dataset[i].z << " " << dataset[i].x << " " << dataset[i].y << " " << dataset[i].cluster << "\n";				
				
				m_labPoint.push_back(dataset[i]);
			}
			else if (dataset[i].pointType == 1)//noise
			{
				Fclustering << dataset[i].id << " " << dataset[i].z << " " << dataset[i].x << " " << dataset[i].y << " " << dataset[i].cluster << "\n";
				m_labPoint.push_back(dataset[i]);
			}
		}

		for (int i = 0; i < corePoint.size(); i++)
		{
			Fclustering << corePoint[i].id << " " << corePoint[i].z << " " << corePoint[i].x << " " << corePoint[i].y << " " << corePoint[i].cluster << "\n";
			m_labPoint.push_back(corePoint[i]);
		}
		Fclustering.close();
	}
	dataset = m_labPoint;
	corePoint.clear();
}
