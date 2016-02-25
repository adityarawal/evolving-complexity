#include "fraserMI.h"

void read_file_data(char *fname, std::vector <double> &X) {

	std::ifstream XtrainFile(fname,std::ios::in);
	double d;
	std::string lineData;

	while (getline(XtrainFile, lineData)) {
			std::stringstream lineStream(lineData);

			while ((lineStream >> d)) {//Read only as many features as there are input nodes
					X.push_back(d);
			}
	}
	XtrainFile.close();

}

/*COMPUTE ENTROPY*/
double fraser_entropy(const vector<double> &X, int feature_index){ //Ranges between 0-max_entropy

	int n;    //Number of points in the feature vector
	int  *Ix; //Index
	double H; //Entropy
	
        //read_file_data("/tmp/x.txt", X); //(X.size() = n+1. Feature vector already consists of unused feature number at the start)
        
	n = X.size();
	Ix = (int *) malloc(sizeof(int) * (n+1));
	quantifyRepeat(&X[0], Ix, n, feature_index);
	H = entropyFraser(Ix, n);
        
        return H;
}


/*COMPUTE MUTUAL INFORMATION*/
double fraser_mutual_information(const vector<double> &X, const vector<double> &Y){

	int n;          //Number of points in the feature vector
	int  *Ix, *Iy;
	int max;
	int *f;
        double minfoFF; //Mutual Information between X and Y

        //read_file_data("/tmp/x.txt", X); //(X.size() = n+1. Feature vector already consists of unused feature number at the start)
	//read_file_data("/tmp/y.txt", Y);

        n = X.size();
	Ix = (int *) malloc(sizeof(int) * (n+1));
	quantifyRepeat(&X[0], Ix, n, 1);

	Iy = (int *) malloc(sizeof(int) * (n+1));
	quantifyRepeat(&Y[0], Iy, n, 2);
	
        f = (int *) malloc(sizeof(int) * (n+1));
	connFeat(Ix, Iy, f, n);
	max = last(f, n);
   // for (int i=1; i<=n;i++){
   // 	std::cout<<X[i-1]<<" "<<Y[i-1]<<" "<<Ix[i]<<" "<<Iy[i]<<" "<<f[i]<<std::endl;
   // }
	minfoFF = mutualInfoFF(f, 1, n, 0, max);

        return minfoFF;
}
