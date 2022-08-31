#include <iostream>
#include <fstream>
#include<cmath>
using namespace std;

int dual_gradientx(int ***rgb, int C,int i,int j,int H,int W)
{
    int val;
    for(int k=0;k<C;k++)
    {
        if(i==0)
        val=(rgb[H-1][j][k]-rgb[(i+1)%H][j][k])*(rgb[H-1][j][k]-rgb[(i+1)%H][j][k]);
        else
        val=(rgb[i-1][j][k]-rgb[(i+1)%H][j][k])*(rgb[i-1][j][k]-rgb[(i+1)%H][j][k]);
    }
    return val;
}
int dual_gradienty(int ***rgb, int C,int i,int j,int H,int W)
{
    int val;
    for(int k=0;k<C;k++)
    {
        if(j==0)
        val=(rgb[i][W-1][k]-rgb[i][(j+1)%W][k])*(rgb[i][W-1][k]-rgb[i][(j+1)%W][k]);
        else
        val=(rgb[i][j-1][k]-rgb[i][(j+1)%W][k])*(rgb[i][j-1][k]-rgb[i][(j+1)%W][k]);
    }
    return val;
}
float ** cal_energy(int ***rgb, int H, int W, int C, int H_, int W_, int C_)
{
    float **energy;
    energy = new float *[H];
    for(int i=0;i<H;i++)
    {
        energy[i] = new float [W];
        for(int j=0;j<W;j++)
        {
            energy[i][j]=sqrt(dual_gradientx(rgb,C,i,j,H,W)+dual_gradienty(rgb,C,i,j,H,W));
        }
    }
return energy;
}
void cal_minPathH(float **energy, int H, int W)
{
    for(int i=1;i<H;i++)
    {
        for(int j=0;j<W;j++)
        {
            if(j==0)
            {
                energy[i][j]=min(energy[i-1][j],energy[i-1][j+1])+energy[i][j];
            }
            else if(j==W-1)
            {
                energy[i][j]=min(energy[i-1][j-1],energy[i-1][j])+energy[i][j];
            }
            else
            {
                energy[i][j]=min(min(energy[i-1][j-1],energy[i-1][j]),energy[i-1][j+1])+energy[i][j];
            }
        }
    }
}
void cal_minPathW(float **energy, int H, int W)
{
    for(int i=0;i<H;i++)
    {
        for(int j=1;j<W;j++)
        {
            if(i==0)
            {
                energy[i][j]=min(energy[i][j-1],energy[i+1][j-1])+energy[i][j];
            }
            else if(i==H-1)
            {
                energy[i][j]=min(energy[i-1][j-1],energy[i][j-1])+energy[i][j];
            }
            else
            {
                energy[i][j]=min(min(energy[i-1][j-1],energy[i][j-1]),energy[i+1][j-1])+energy[i][j];
            }
        }
    }
}
int* path_(float **energy, int H, int W)
{
    int * path_idx=new int [H];
    float min = 2147483647;
        for(int j=W-1;j>=0;j--)
        {
            if(min>energy[H-1][j])
            {
                min=energy[H-1][j];
                path_idx[H-1]=j;
            }
        }
        for(int i=H-2;i>=0;i--)
        {
            int j=path_idx[i+1];
            if(j==0)
            {
                if(energy[i][j]<energy[i][j+1])
                path_idx[i]=j;
                else
                path_idx[i]=j+1;
            }
            else if(j==W-1)
            {
                if(energy[i][j]<energy[i][j-1])
                path_idx[i]=j;
                else
                path_idx[i]=j-1;
            }
            else
            {
                if(energy[i][j]<energy[i][j+1])
                {
                    if(energy[i][j]<energy[i][j-1])
                    path_idx[i]=j;
                    else
                    path_idx[i]=j-1;
                }
                else
                {
                    if(energy[i][j-1]<energy[i][j+1])
                    path_idx[i]=j-1;
                    else
                    path_idx[i]=j+1;
                }
            }
        }
        return path_idx;
}
int* path1_(float **energy, int H, int W)
{
    int * path_idx=new int [W];
    float min = 2147483647;
        for(int j=H-1;j>=0;j--)
        {
            if(min>energy[j][W-1])
            {
                min=energy[j][W-1];
                path_idx[W-1]=j;
            }
        }
        for(int j=W-2;j>=0;j--)
        {
            int i=path_idx[j+1];
            if(i==0)
            {
                if(energy[i][j]<energy[i+1][j])
                path_idx[j]=i;
                else
                path_idx[j]=i+1;
            }
            else if(i==H-1)
            {
                if(energy[i][j]<energy[i-1][j])
                path_idx[j]=i;
                else
                path_idx[j]=i-1;
            }
            else
            {
                if(energy[i][j]<energy[i+1][j])
                {
                    if(energy[i][j]<energy[i-1][j])
                    path_idx[j]=i;
                    else
                    path_idx[j]=i-1;
                }
                else
                {
                    if(energy[i-1][j]<energy[i+1][j])
                    path_idx[j]=i-1;
                    else
                    path_idx[j]=i+1;
                }
            }
        }
        return path_idx;
}
void reduceW_rgb(int ***rgb, int* path_idx, int H, int W, int C, int H_, int W_, int C_)
{
    for(int i=0;i<H;i++)
    {
        for(int j=0;j<W-1;j++)
        {
            if(j>=path_idx[i])
            {
                for(int k=0;k<C;k++)
                {
                    rgb[i][j][k]=rgb[i][j+1][k];
                }
            }
        }
    }
}
void reduceH_rgb(int ***rgb, int* path_idx, int H, int W, int C, int H_, int W_, int C_)
{
    for(int j=0;j<W;j++)
    {
        for(int i=0;i<H-1;i++)
        {
            if(i>=path_idx[j])
            {
                for(int k=0;k<C;k++)
                {
                    rgb[i][j][k]=rgb[i+1][j][k];
                }
            }
        }
        
    }
}
void solve(int ***rgb, int H, int W, int C, int H_, int W_, int C_) {
    // We've provided you the driver.py and main.cpp for your convinience
    // Please go through them and understand the file handling and input/output format
    // Feel free to experiment on your own

    // can access the array using rgb[i][j][k] like a regular 3D array

    // Write your code here
    float **energy;
    int *path_idx;
    for(int p=0;p<W-W_;p++){
    energy = cal_energy(rgb, H, W-p, C, H_, W_, C_);
    cal_minPathH(energy,H,W-p);
    path_idx=path_(energy,H,W-p);
    reduceW_rgb(rgb,path_idx, H, W-p, C, H_, W_, C_);
    }
    W=W_;
    for(int p=0;p<H-H_;p++)
    {
        energy = cal_energy(rgb, H-p, W, C, H_, W_, C_);
        cal_minPathW(energy,H-p,W);
        path_idx=path1_(energy,H-p,W);
        reduceH_rgb(rgb,path_idx, H-p, W, C, H_, W_, C_);
    }
}
int main() {
    string ip_dir = "/home/dell/Downloads/DSAPS/Q3/data/input/";
    string ip_file = "rgb_in.txt";
    ifstream fin(ip_dir + ip_file);

    int H, W, C;
    fin >> H >> W >> C;

    int ***rgb;
    rgb = new int **[H];
    for(int i = 0; i < H; ++i) {
        rgb[i] = new int *[W];
        for(int j = 0; j < W; ++j) {
            rgb[i][j] = new int[C];
            for(int k = 0; k < C; ++k)
                fin >> rgb[i][j][k];
        }
    }
    fin.close();

    int H_, W_, C_;
    cout << "Enter new value for H (must be less than " << H << "): ";
    cin >> H_;
    cout << "Enter new value for W (must be less than " << W << "): ";
    cin >> W_;
    cout << '\n';
    C_ = C;

    solve(rgb, H, W, C, H_, W_, C_);

    string op_dir = "/home/dell/Downloads/DSAPS/Q3/data/output/";
    string op_file = "rgb_out.txt";
    ofstream fout(op_dir + op_file);
    
    fout << H_ << " " << W_ << " " << C_ << endl;
    for(int i = 0; i < H_; ++i) {
        for(int j = 0; j < W_; ++j) {
            for(int k = 0; k < C_; ++k) {
                fout << rgb[i][j][k] << " ";
            }
        }
        fout << '\n';
    }
    fout.close();
    return 0;
}