struct metacolor
{
	int colorID;
	float H;
	float S;
	float V;
	float R;
	float G;
	float B;
	float L;
	float A;
	float Bright;
	float weight;
	int clusterID;
	bool selected = false;
};
struct BeeGroup
{
	std::vector<metacolor> code;  
	double trueFit;//Record the minimum  
	double fitness;
	double rfitness;//fitness ratio 
	int trail;//number of experiments  
}*Bee;
BeeGroup BestSource;//Record the best source 
std::vector<metacolor> inputColor;

const int NP = 40;//Population size  
const int FoodNumber = NP / 2;//Bee size  
const int limit = 20;//Beyond this limit, the bee has not been renewed to become a scout bee 
const int maxCycle = 1000;  

/*****parameters*****/
int D = 11;// number of color  
const double lb = -100;//Lower bound  
const double ub = 100;//upper bound
const float lbH = 0;//Lower bound of H(HSV)    
const float ubH = 360;//upper bound of H(HSV) 
const float lbV = 25;//Lower bound of V(HSV)    
const float lbS = 10;//Lower bound of S(HSV)     
const float ubSV = 100;//upper bound of SV(HSV)
const float lbL = 0;//Lower bound of L(lab)     
const float ubL = 100;//upper bound of L(lab) 
const float lbAB = -128;//Lower bound of ab(lab)  
const float ubAB = 127;//upper bound of ab(lab) 
double result[10000] ;
int limitDis = 25;  //Minimum distance between two colors
int limitDisAss = 50;  //Maximum distance of related colors
int limitDistoIni = 60;  //Minimum distance between ini-colors
float DBSCANEps = 30;	

//random function
double  random(double start, double end)
{
	return start + (end - start)*rand() / (RAND_MAX + 1.0);
}

void TrainABC()
{
	initilizeBee();
	initilize();  
	
	MemorizeBestSource();//save best source  

	int gen = 0; 
	bool first = true; bool first1 = true; bool first2 = true;
	while (gen<maxCycle)
	{
		sendEmployedBees();

		CalculateProbabilities();

		sendOnlookerBees();

		MemorizeBestSource();

		sendScoutBees();

		MemorizeBestSource();

		gen++;
	}
}

void initilize() 
{
	int i, j;
	
	for (i = 0; i<FoodNumber; i++)
	{
		NectarSource[i].code.resize(D);
		OnLooker[i].code.resize(D);
		EmployedBee[i].code.resize(D);

		for (j = 0; j<D; j++)
		{
			NectarSource[i].code[j] = inputColor[j];
			metacolor color = NectarSource[i].code[j];
			int clustertype = colorCluster[color.clusterID].type;
			while (true)
			{
				if (clustertype != 1)
				{
					NectarSource[i].code[j].H = inputColor[j].H + random(-10, 10);
				}
				NectarSource[i].code[j].S = inputColor[j].S + random(-20, 20);
				NectarSource[i].code[j].V = inputColor[j].V + random(-15, 15);
				colorAdjustHSV(NectarSource[i].code[j]);
				calmetaColorHSV(NectarSource[i].code[j]);				
				if (calcolorDis(NectarSource[i].code[j], inputColor[j]) < limitDistoIni)
				{
					break;
				}			
			}
			
			EmployedBee[i].code[j] = NectarSource[i].code[j];
			OnLooker[i].code[j] = NectarSource[i].code[j];
			BestSource.code[j] = NectarSource[0].code[j];
		}
		/****Initialization of the nectar source*****/
		NectarSource[i].trueFit = calculationTruefit(NectarSource[i]);
		NectarSource[i].fitness = calculationFitness(NectarSource[i].trueFit);
		NectarSource[i].rfitness = 0;
		NectarSource[i].trail = 0;
		/****Initialization of the employed bees*****/
		EmployedBee[i].trueFit = NectarSource[i].trueFit;
		EmployedBee[i].fitness = NectarSource[i].fitness;
		EmployedBee[i].rfitness = NectarSource[i].rfitness;
		EmployedBee[i].trail = NectarSource[i].trail;
		/****Initialization of the onlooker bees****/
		OnLooker[i].trueFit = NectarSource[i].trueFit;
		OnLooker[i].fitness = NectarSource[i].fitness;
		OnLooker[i].rfitness = NectarSource[i].rfitness;
		OnLooker[i].trail = NectarSource[i].trail;
	}
	/*****Initialization of the best  food source*****/
	BestSource.trueFit = NectarSource[0].trueFit;
	BestSource.fitness = NectarSource[0].fitness;
	BestSource.rfitness = NectarSource[0].rfitness;
	BestSource.trail = NectarSource[0].trail;

	
}

void initilizeBee()
{
	if (Bee)
	{
		delete[] Bee;
		Bee = 0;
	}
	Bee = new BeeGroup[FoodNumber];
	NectarSource=new BeeGroup[FoodNumber];
	EmployedBee = new BeeGroup[FoodNumber]; 
	OnLooker = new BeeGroup[FoodNumber];

	BestSource.rfitness = 0;
	BestSource.fitness = 0;
	BestSource.code.resize(D);
	BestSource.trueFit = 0;
	BestSource.trail = 0;

	iniEnergy = 0.0f;
	for (size_t i = 0; i < D; i++)
	{
		iniEnergy += inputColor[i].weight*(inputColor[i].R + inputColor[i].G + inputColor[i].B);
	}
}

//Calculate hsv&lab by input map rgb
void calmetaColor(metacolor &c)
{
	vector<double> lab =colorTransfer::RGB2Lab(c.R,c.G,c.B);
	vector<double> hsv= colorTransfer::RGB2HSV(c.R, c.G, c.B);

	c.L = lab[0];
	c.A = lab[1];
	c.Bright = lab[2];

	c.H= hsv[0];
	c.S = hsv[1];
	c.V = hsv[2];
}

//Judging hsv color range
void colorAdjustHSV(metacolor &c)
{
	if (c.H > ubH)
	{
		c.H = c.H - ubH;
	}
	if (c.H < lbH)
	{
		c.H = c.H + ubH;
	}
	if (c.S > ubSV)
	{
		c.S = ubSV;
	}
	if (c.S < lbS)
	{
		c.S = lbS;
	}
	if (c.V > ubSV)
	{
		c.V = ubSV;
	}
	if (c.V < lbV)
	{
		c.V = lbV;
	}
}

//Judging lab color range
void colorAdjustLab(metacolor &c)
{
	if (c.L > ubL)
	{
		c.L = ubL;
	}
	if (c.L < lbL)
	{
		c.L = lbL;
	}
	if (c.A > ubAB)
	{
		c.A = ubAB;
	}
	if (c.A < lbAB)
	{
		c.A = lbAB;
	}
	if (c.Bright > ubAB)
	{
		c.Bright = ubAB;
	}
	if (c.Bright < lbAB)
	{
		c.Bright = lbAB;
	}
}

//Calculate rgb&lab by input hsv
void calmetaColorHSV(metacolor &c)
{
	vector<double> lab = colorTransfer::HSB2Lab(c.H, c.S, c.V);


	double hsb[] = { c.H, c.S, c.V };
	int rgb[3];
	colorTransfer::HSB2RGB(hsb,rgb);

	c.L = lab[0];
	c.A = lab[1];
	if (lab[2] == 127)
	{
		int iii = 0;
	}
	c.Bright = lab[2];

	c.R = rgb[0];
	c.G = rgb[1];
	c.B = rgb[2];

	if (c.R == 0 &&c.G==0)
		int xx = 111;
}

//Calculate rgb&hsv by input lab
void calmetaColorLAB(metacolor &c)
{
	vector<double> hsv = colorTransfer::lab2hsv(c.L, c.A, c.Bright);
	vector<double> rgb = colorTransfer::LabToRgb(c.L, c.A, c.Bright);

	c.H= hsv[0];
	c.S = hsv[1];
	c.V = hsv[2];

	c.R = rgb[0];
	c.G = rgb[1];
	c.B = rgb[2];

	if (c.R == 0 && c.G == 0)
		int xx = 111;
}

//calculate Truefit
double calculationTruefit(BeeGroup bee)  
{
	double truefit = 0;
		
	float Energy = 0.0;

	for (int i = 0; i < bee.code.size(); i++)
	{
		Energy += (bee.code[i].G + bee.code[i].R + bee.code[i].B)*bee.code[i].weight;
	}


	if (Energy < 0)
		int xx = 11;
	truefit = Energy / iniEnergy;

	return truefit;
}

//calculate Fitness
double calculationFitness(double truefit)  
{
	double fitnessResult = 0;
	if (truefit >= 0)
	{
		fitnessResult = 1 / (truefit + 1);
	}
	else
	{
		fitnessResult = 1 + abs(truefit);
	}
	return fitnessResult;
}

//calculate color distance
float calcolorDis(const metacolor &c1, const metacolor &c2)
{
	float distant = sqrt(pow((c1.L - c2.L), 2) + pow((c1.A - c2.A), 2) + pow((c1.Bright - c2.Bright), 2));
	return distant;
}

//Judging color semantic relations
bool judegeColorRelu(int id,BeeGroup &bee1)
{
	metacolor color = bee1.code[id];
	int clustertype = colorCluster[color.clusterID].type; 	//0 is differentiation,1 is order,2 is association

	if (calcolorDis(bee1.code[id], inputColor[id])>limitDistoIni)
		return false;

	for (int i = 0; i < D; i++)
	{
		if (i == id)
			continue;
		if(calcolorDis(bee1.code[i],color)<limitDis)
			return false;
	}

	if (clustertype == 0)
	{
		vector<int> cs = colorCluster[color.clusterID].colorIDs;
		for (int i = 0; i < cs.size(); i++)
		{
			if (cs[i] != id)
			{
				if (calcolorDis(bee1.code[i], color)>limitDis)
					return false;
			}
		}
	}
	else if (clustertype == 1)
	{
		vector<int> cs = colorCluster[color.clusterID].colorIDs;
		for (int i = 0; i < cs.size(); i++)
		{
			if (cs[i] != id)
			{
				float iii = calcolorDis(bee1.code[i], inputColor[i]);
				if (calcolorDis(bee1.code[i], inputColor[i])>limitDis)
					return false;
			}
		}
	}
	else if (clustertype == 2)
	{
		vector<int> cs = colorCluster[color.clusterID].colorIDs;
		for (int i = 0; i < cs.size(); i++)
		{
			if (cs[i] != id)
			{
				if (calcolorDis(bee1.code[i], color)>limitDisAss)
					return false;
			}
		}
	}
}

bool judegeColorSelected(int id, BeeGroup &bee)
{
	if (colorCluster[bee.code[id].clusterID].type == 1 && !bee.code[id].selected)
	{
		return false;
	}
	return true;
}

vector<float> getInterParameterLine(double y1, double y2, int count)
{
	float a, b;

	float x1 = 1;
	float x2 = count;

	a = (y2 - y1) / (x2 - x1);
	b = y1 - a * x1;

	return vector<float>{a, b};
}

void interpColor(int id, BeeGroup &bee)
{
	vector<metacolor> orcolors;
	int clusterid = bee.code[id].clusterID;
	int colornum = colorCluster[clusterid].colorIDs.size();
	
	for (int i = 0; i < colornum; i++)
	{
		orcolors.emplace_back(bee.code[colorCluster[clusterid].colorIDs[i]]);
	}
	std::sort(orcolors.begin(), orcolors.end(), [](const metacolor& x, const metacolor& y)
	{
		return x.L > y.L;
	});
	vector<metacolor> orcolors1;
	for(int i = 0; i < colornum; i++)
	{
		if (bee.code[orcolors[i].colorID].selected)
			orcolors1.emplace_back(orcolors[i]);
	}
	
	vector<float> parameterS = getInterParameterLine(orcolors1[0].S, orcolors1[2].S, colornum);
	vector<float> parameterV = getInterParameterLine(orcolors1[0].V, orcolors1[2].V, colornum);

	std::vector<double> tt;
	for (int j = 1; j < colornum+1 ; j++)
	{
		
		int count = 50;
		while (true)
		{
			double Rij = random(0, 1);
			bee.code[orcolors[j - 1].colorID].S = parameterS[0] * j + parameterS[1];
			bee.code[orcolors[j - 1].colorID].V = parameterV[0] * j + parameterV[1];

			colorAdjustHSV(bee.code[orcolors[j - 1].colorID]);
			calmetaColorHSV(bee.code[orcolors[j - 1].colorID]);
			if ((count--)< 0)
			{
				break;
			}
			for (int t = 0; t < D; t++)
			{
				if (t == (j-1))
					continue;
				if (calcolorDis(bee.code[orcolors[t].colorID], bee.code[orcolors[j - 1].colorID]) > limitDis)
					break;
			}
			if (calcolorDis(bee.code[orcolors[j - 1].colorID], inputColor[j - 1]) < limitDistoIni)
			{
				break;
			}
		}
	}

}

void ChangeBee(int id,BeeGroup &bee1,const BeeGroup &bee2, const BeeGroup &bee3)
{
	metacolor inicolor= bee1.code[id];

	double Rij = random(-1, 1);
	
	if (!judegeColorSelected(id, bee1))
		return;

	if (colorCluster[inicolor.clusterID].type != 1)
	{
		bee1.code[id].H = bee2.code[id].H + Rij*(bee2.code[id].H - bee3.code[id].H);
		if (bee1.code[id].H >ubH)
		{
			bee1.code[id].H = bee1.code[id].H - ubH;
		}
		if (bee1.code[id].H <lbH)
		{
			bee1.code[id].H = bee1.code[id].H + ubH;
		}
	}
	Rij = random(-1, 1);
	bee1.code[id].S = bee2.code[id].S + Rij*(bee2.code[id].S - bee3.code[id].S);
	if (bee1.code[id].S >ubSV)
	{
		bee1.code[id].S = ubSV;
	}
	if (bee1.code[id].S <lbS)
	{
		bee1.code[id].S = lbS;
	}
	Rij = random(-1, 1);
	bee1.code[id].V = bee2.code[id].V + Rij*(bee2.code[id].V - bee3.code[id].V);
	if (bee1.code[id].V >ubSV)
	{
		bee1.code[id].V = ubSV;
	}
	if (bee1.code[id].V <lbV)
	{
		bee1.code[id].V = lbV;
	}

	calmetaColorHSV(bee1.code[id]);

	if (calcolorDis(bee1.code[id], inputColor[id]) > limitDistoIni)
	{
		bee1.code[id] = inicolor;
		return;
	}

	if (colorCluster[inicolor.clusterID].type == 1)
	{
		interpColor(id, bee1);
		for (int i = 0; i < colorCluster[inicolor.clusterID].colorIDs.size() - 1; i++)
			for (int j = i; j < colorCluster[inicolor.clusterID].colorIDs.size(); j++)
			{
				if (calcolorDis(bee1.code[colorCluster[inicolor.clusterID].colorIDs[i]], bee1.code[colorCluster[inicolor.clusterID].colorIDs[j]])<limitDis)
				{
					bee1.code[id] = inicolor;
					return;
				}
			}
	}

	if (!judegeColorRelu(id, bee1))
	{
		bee1.code[id] = inicolor;
		return;
	}
}

void sendEmployedBees()  
{
	int i, j, k;
	int param2change; 
	double Rij = random(-1, 1);//random 
	for (i = 0; i<FoodNumber; i++)
	{

		param2change = (int)random(0, D);  
		/******i!=k********/
		while (1)
		{
			k = (int)random(0, FoodNumber);
			if (k != i)
			{
				break;
			}
		}
		for (j = 0; j<D; j++)
		{
			EmployedBee[i].code[j] = NectarSource[i].code[j];
		}

		ChangeBee(param2change, EmployedBee[i], NectarSource[i], NectarSource[k]);
	
		EmployedBee[i].trueFit = calculationTruefit(EmployedBee[i]);
		EmployedBee[i].fitness = calculationFitness(EmployedBee[i].trueFit);

		/******Greedy selection strategy*******/
		if (EmployedBee[i].trueFit<NectarSource[i].trueFit)
		{
			for (j = 0; j<D; j++)
			{
				NectarSource[i].code[j] = EmployedBee[i].code[j];
			}
			NectarSource[i].trail = 0;
			NectarSource[i].trueFit = EmployedBee[i].trueFit;
			NectarSource[i].fitness = EmployedBee[i].fitness;
		}
		else
		{
			NectarSource[i].trail++;
		}
	}
}

void CalculateProbabilities()  
{
	int i;
	double maxfit;
	maxfit = NectarSource[0].fitness;
	for (i = 1; i<FoodNumber; i++)
	{
		if (NectarSource[i].fitness>maxfit)
			maxfit = NectarSource[i].fitness;
	}

	for (i = 0; i<FoodNumber; i++)
	{
		NectarSource[i].rfitness = (0.9*(NectarSource[i].fitness / maxfit)) + 0.1;
	}
}

void sendOnlookerBees()  
{
	int i, j, t, k;
	double R_choosed;  
	int param2change;  
	double Rij; 
	i = 0;
	t = 0;
	while (t<FoodNumber)
	{

		R_choosed = random(0, 1);
		if (R_choosed<NectarSource[i].rfitness)//according to the probability of being selected  
		{
			t++;
			param2change = (int)random(0, D);

			/******i!=k********/
			while (1)
			{
				k = (int)random(0, FoodNumber);
				if (k != i)
				{
					break;
				}
			}
			for (j = 0; j<D; j++)
			{
				OnLooker[i].code[j] = NectarSource[i].code[j];
			}

			ChangeBee(param2change, OnLooker[i], NectarSource[i], NectarSource[k]);

			OnLooker[i].trueFit = calculationTruefit(OnLooker[i]);
			OnLooker[i].fitness = calculationFitness(OnLooker[i].trueFit);

			/****Greedy selection strategy******/
			if (OnLooker[i].trueFit<NectarSource[i].trueFit)
			{
				for (j = 0; j<D; j++)
				{
					NectarSource[i].code[j] = OnLooker[i].code[j];
				}
				NectarSource[i].trail = 0;
				NectarSource[i].trueFit = OnLooker[i].trueFit;
				NectarSource[i].fitness = OnLooker[i].fitness;
			}
			else
			{
				NectarSource[i].trail++;
			}
		}
		i++;
		if (i == FoodNumber)
		{
			i = 0;
		}
	}
}

/*******only one scout bee**********/
void sendScoutBees()
{
	int maxtrialindex, i, j;
	double R; 
	maxtrialindex = 0;
	for (i = 1; i<FoodNumber; i++)
	{
		if (NectarSource[i].trail>NectarSource[maxtrialindex].trail)
		{
			maxtrialindex = i;
		}
	}
	if (NectarSource[maxtrialindex].trail >= limit)
	{
		/*******reinitialize*********/
		for (j = 0; j<D; j++)
		{
			if (!judegeColorSelected(j, NectarSource[maxtrialindex]))
				continue;

			auto inicolor = NectarSource[maxtrialindex].code[j];

			if (colorCluster[inicolor.clusterID].type != 1)
			{
				R = random(0, 1);
				NectarSource[maxtrialindex].code[j].H = lbH + R*(ubH - lbH);
			}
			R = random(0, 1);
			NectarSource[maxtrialindex].code[j].S = lbS + R*(ubSV - lbS);
			R = random(0, 1);
			NectarSource[maxtrialindex].code[j].V = lbV + R*(ubSV - lbV);

			calmetaColorHSV(NectarSource[maxtrialindex].code[j]);

			if (calcolorDis(NectarSource[maxtrialindex].code[j], inputColor[j]) > limitDistoIni)
			{
				NectarSource[maxtrialindex].code[j] = inicolor;
				continue;
			}

			if (!judegeColorRelu(j, NectarSource[maxtrialindex]))
			{
				NectarSource[maxtrialindex].code[j] = inicolor;
				continue;
			}

			if (colorCluster[inicolor.clusterID].type == 1)
			{
				interpColor(j, NectarSource[maxtrialindex]);
				for (int i0 = 0; i0 < colorCluster[inicolor.clusterID].colorIDs.size() - 1; i0++)
					for (int j0 = i0; j0 < colorCluster[inicolor.clusterID].colorIDs.size(); j0++)
					{
						if (calcolorDis(NectarSource[maxtrialindex].code[colorCluster[inicolor.clusterID].colorIDs[i0]], 
							NectarSource[maxtrialindex].code[colorCluster[inicolor.clusterID].colorIDs[j0]]) < limitDis)
						{
							NectarSource[maxtrialindex].code[j] = inicolor;
							continue;
						}
					}
			}
			
		}
		NectarSource[maxtrialindex].trail = 0;
		NectarSource[maxtrialindex].trueFit = calculationTruefit(NectarSource[maxtrialindex]);
		NectarSource[maxtrialindex].fitness = calculationFitness(NectarSource[maxtrialindex].trueFit);
	}
}

//Memory the best sources
void MemorizeBestSource() 
{
	int i, j;
	for (i = 1; i<FoodNumber; i++)
	{
		if (NectarSource[i].trueFit<BestSource.trueFit)
		{
			for (j = 0; j<D; j++)
			{
				BestSource.code[j] = NectarSource[i].code[j];
			}
			BestSource.trueFit = NectarSource[i].trueFit;
		}
	}
}
