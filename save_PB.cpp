#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm> 
#include <cmath>
using namespace std;

#ifndef BYTE
typedef unsigned char BYTE;
#endif
#ifndef WORD
typedef unsigned short WORD;
#endif
#ifndef DWORD
typedef unsigned long DWORD;
#endif

//#define MAXN 2048	// You can increase this value if N>1024.
#define MAXN 1024	// You can increase this value if N>1024.

#define W 32
#define MAXDWSZ ((MAXN/W) + 1)

#define BLK_SZ 8
#define MAXMBLK_NUM	(MAXN/BLK_SZ+2)
#define N 32

DWORD f[MAXDWSZ] = {0}; 

int DWSZ, BLK_NUM, Highest_time_of_f;

BYTE TBL[65536], TBH[65536];

int Nowf[30]={0x11B,0x11D,0x12B,0x12D,0x139,0x13F,0x14D,0x15F,0x163,0x165,0x169,0x171,0x177,0x17B,0x187,0x18B,0x18D,0x19F,0x1A3,0x1A9,0x1B1,0x1BD,0x1C3,0x1CF,0x1D7,0x1DD,0x1E7,0x1F3,0x1F5,0x1F9};


#define SET_BASIC_PARAMETERS	{ 	f[0] = 0x1c3 ; Highest_time_of_f = 8;					\
								if ( 0 != (N%W) ) DWSZ = ((N/W) + 1); else DWSZ=N/W;	\
								if ( 0 == (N % BLK_SZ) )								\
								  BLK_NUM = ((N/BLK_SZ));								\
								  else BLK_NUM = ((N/BLK_SZ)+1);						\
								}

#define ln_copy(a, b)  for (int zqa = 0; zqa < DWSZ; zqa++) a[zqa] = b[zqa];

void precompute_table()
{
	int i;
	DWORD bh, bl;
	DWORD w, wt;

	for(bl = 0; bl < 256; bl++) 
		for(bh = bl; bh < 256; bh++) 
		{
			bl &= 0xff; bh &= 0xff;
			w = 0; wt = (DWORD)bl;
			for (i = 0; i < 8; i++) if (1 & (bh>>i) ) w ^= (wt<<i);;	// w=bl*bh

			wt = ((wt << 8) | ((DWORD)bh)) & 0xffff;
			TBL[wt] = ((BYTE)(w & 0xff)); 
			TBH[wt] = ((BYTE)(w >> 8)); 
			wt = ((wt << 8) | ((DWORD)bl)) & 0xffff;
			TBL[wt] = ((BYTE)(w & 0xff)); 
			TBH[wt] = ((BYTE)(w >> 8)); 
		}
}

void poly_mul(DWORD *c, DWORD *a, DWORD *b)
{
	BYTE r[4*2*MAXDWSZ+1];
	int zqa, i, j, h;
	BYTE *ba = (BYTE*)a, *bb = (BYTE*)b;
	DWORD *pd, k;
	pd = (DWORD*)r;
	for(i = 0; i < 2 * MAXDWSZ; i++) pd[i] = 0;

	for(i = 0; i < BLK_NUM; i++)
		for(j = 0; j < BLK_NUM; j++){
			k = ((DWORD)ba[i]) | (((DWORD)bb[j]) << 8); 
			r[i+j] ^= TBL[k]; r[i+j+1] ^= TBH[k];
		}
	pd = (DWORD*)r;
	ln_copy(c, pd);
}

void poly_mod(DWORD *c, DWORD *ff)
{
	int shift_bits = 0;
	for(int i = BLK_NUM - 1; i >= 0; i--){
		for(int j = 7; j >= 0; j--){
			BYTE *bc = (BYTE*)c;
			if((((DWORD)bc[i]) >> j) & 1){
				shift_bits = j + 8 * i - Highest_time_of_f;
				if(shift_bits < 0) return;
				(*c) ^= ((*ff) << shift_bits);
			}
		}
	}	
}

#define Matrix_Msize 10
#define Matrix_size 4
#define Reverse_Mtrx_Hmwt 0
#define txt_name "testc_1c3.txt" 

//87610 
int mul[260][260],exprAB[260][8][8][8],now[260][8][8],ny[260];


void genexpr(int R)
{
	int mod=f[0];
    for(int i=0;i<256;i++)
        for(int j=0;j<256;j++)
        {
            int x=i,y=j,z=0;
            while(y)
            {
                if(y&1)
                    z^=x;
                x=(x<<1)^((((x<<1)&256))?mod:0);
                y>>=1;
            }
            mul[i][j]=z;
            if(z==1)
                ny[i]=j;
        }
    for(int i=0;i<8;i++)
        for(int j=0;j<8;j++)
            for(int k=0;k<8;k++)
                if(mul[ny[R]][mul[1<<i][1<<j]]&(1<<k))
                {
                    exprAB[R][k][j][now[R][k][j]]=i;
                    now[R][k][j]++;
                }
    for(int i=0;i<8;i++)
        for(int j=0;j<8;j++)
            exprAB[R][i][j][now[R][i][j]]=-1;
}
int coefficient[4][8] = {0};

void Initialize_c()
{
	for(int i = 0; i < 4; i++)
		for(int j = 0; j < 8; j++)
		{
			coefficient[i][j] = 0;
		}
}

struct aa
{
    int hmw, rhmw;
    int XOR_num_in, XOR_num_out;
    int XOR_delay_in, XOR_delay_out;
    int XOR_num_all, XOR_delay_all;
    string vec_in, vec_out;
}pb_vec[260][32768];

bool cmp(const aa &x, const aa &y) 
{
	if(x.XOR_delay_all == y.XOR_delay_all) 
    	return x.XOR_num_all < y.XOR_num_all; //从小到大排<，若要从大到小排则> 
    else return x.XOR_delay_all < y.XOR_delay_all;
} 

int cnt = 0;

int count_XOR_delay()
{
	int max_delay = 0, tmp_delay = 0;
	for(int j = 0; j < 8; j++)
	{
		tmp_delay = 0;
		for(int i = 0; i < 4; i++)
		{
			tmp_delay += coefficient[i][j];
		}
		// cout << tmp_delay << endl;
		if(max_delay < tmp_delay) max_delay = tmp_delay;
	}
	double tmp = log(max_delay)/log(2);
	max_delay = ceil(tmp);
	return max_delay;
}

void count_XOR_num(DWORD dw, int order, int is_vec_in, int row, int R)
{
	int arr[8] = {0}; int aaa = dw;
	int w = 0;
	while(dw > 0){
		if(dw & 1) arr[w] = 1;
		dw >>= 1;
		w++;
	}
	int expr_size = 8, ele_num = 8, XOR_num = 0;
	for(int i = 0; i < expr_size; i++)
	{
		int XOR_num_per_row = 0;
		for(int j = 0; j < expr_size; j++)
		{
			int XOR_num_per_ele = 0;
			for(int k = 0; k < ele_num; k++)
			{
				if(exprAB[R][i][j][k] == -1) break;
				if(arr[exprAB[R][i][j][k]] == 1) XOR_num_per_ele ^= 1;
			}
			XOR_num_per_row += XOR_num_per_ele;
		}
		XOR_num += XOR_num_per_row;
		coefficient[row][i] = XOR_num_per_row;
	}
	XOR_num -= 8;
	if(is_vec_in == 1) pb_vec[R][order].XOR_num_in += XOR_num;
	else if(is_vec_in == 0) pb_vec[R][order].XOR_num_out += XOR_num;	
}

int getA(DWORD arcs[Matrix_Msize][Matrix_Msize], int n)//按第一行展开计算|A|
{
    if(n == 1) return arcs[0][0];
    int ans = 0;
    DWORD temp[Matrix_Msize][Matrix_Msize];
    int i, j, k;
    for(i = 0; i < n; i++)
    {
        for(j = 0; j < n - 1; j++)
        {
            for(k = 0; k < n - 1; k++)
            {
                temp[j][k] = arcs[j + 1][(k >= i) ? k + 1 : k];

            }
        }
        int t = getA(temp, n - 1);
        if(t == 0) return 0;
  //       if(i % 2 == 0)  ans += t * arcs[0][i];
		// else ans -= t * arcs[0][i];
		DWORD c[MAXDWSZ] = {0}; 
		DWORD a[MAXDWSZ] = {0}; a[0] = t;
		DWORD b[MAXDWSZ] = {0}; b[0] = arcs[0][i];
		poly_mul(c, a, b);
		poly_mod(c, f);       
        ans ^= c[0];
    }
    return ans;
}

void getAStart(DWORD arcs[Matrix_Msize][Matrix_Msize], int n, DWORD ans[Matrix_Msize][Matrix_Msize])//计算每一行每一列的每个元素所对应的余子式，组成A*
{
    if(n == 1)
    {
        ans[0][0] = 1;
        return;
    }
    int i, j, k, t;
    DWORD temp[Matrix_Msize][Matrix_Msize];
    for(i = 0; i < n; i++)
    {
        for(j = 0; j < n; j++)
        {
            for(k = 0; k < n - 1; k++)
            {
                for(t = 0; t < n - 1; t++)
                {
                    temp[k][t] = arcs[k >= i ? k + 1 : k][t >= j ? t + 1 : t];
                }
            }

            ans[j][i] = getA(temp, n - 1);
        }
    }
}

int Hamming_weight(DWORD a)
/*return the Hamming_weight of a >0 */
{
	int w = 0;
	while(a > 0){
		if(a & 1) w++;
		a >>= 1;
	}
	return w;
}

void inverse_vec(DWORD input_array[Matrix_size])
{
	int hm_num = 0;
    DWORD arcs[Matrix_Msize][Matrix_Msize];
    for(int i = 0; i < Matrix_size; i++)
    {
    	for(int j = 0; j < Matrix_size; j++)
    	{
    		arcs[i][j] = input_array[(j - i + Matrix_size) % Matrix_size];
    	}
    }	

    DWORD astar[Matrix_Msize][Matrix_Msize];
    int a = getA(arcs, Matrix_size);
    if (a == 0) return;  

	DWORD MI_a[MAXDWSZ] = {0}; DWORD b[MAXDWSZ] = {0}; DWORD output_array[Matrix_size];

    getAStart(arcs, Matrix_size, astar);

    // multiplicative inverse of a
	MI_a[0] = 1;  
	b[0] = a;
    for(int kk = 0; kk < 7; kk++)
    {        	
		poly_mul(b, b, b);   
		poly_mod(b, f);
		poly_mul(MI_a, MI_a, b);   
		poly_mod(MI_a, f);
    }

    //逆矩阵
	for(int i = 0; i < Matrix_size; i++)
	{
		for(int j = 0; j < Matrix_size; j++)
		{
	        b[0] = astar[i][j];
			poly_mul(b, MI_a, b);   
			poly_mod(b, f);
			astar[i][j] = *b;
		}
    }

    for(int j = 0; j < Matrix_size; j++)
    {
        b[0] = astar[0][j];
		hm_num += Hamming_weight(*b);
    }

	for(int R=1;R<=255;R++)
	{
		pb_vec[R][cnt].rhmw = hm_num;	
	    /*   PB_in   */
		MI_a[0] = R; hm_num = 0;
		pb_vec[R][cnt].vec_in += "(";
		Initialize_c();
		for(int i = 0; i < Matrix_size; i++)
		{
			b[0] = input_array[i];	
			poly_mul(b, MI_a, b);   
			poly_mod(b, f);		
			count_XOR_num((*b), cnt, 1, i,R);
		    char str[10];
		    itoa((*b), str, 16);	
		    if(strlen(str) < 2)
		    {
		    	char addzero[10] = "0";
		    	strcat(addzero, str); 
		    	strcpy(str, addzero);
		    }
			pb_vec[R][cnt].vec_in += str;
			hm_num += Hamming_weight(*b);	
			if(i != Matrix_size - 1) pb_vec[R][cnt].vec_in += " ";
		}
		pb_vec[R][cnt].vec_in += ") "; pb_vec[R][cnt].hmw = hm_num;
		pb_vec[R][cnt].XOR_num_in = 4 * (pb_vec[R][cnt].XOR_num_in + 3 * 8);
		pb_vec[R][cnt].XOR_delay_in = count_XOR_delay();

		/*   PB_out   */
		
		pb_vec[R][cnt].vec_out += " (";
	    Initialize_c();
	    for(int j = 0; j < Matrix_size; j++)
	    {
			b[0] = astar[0][j];	
			poly_mul(b, MI_a, b);   
			poly_mod(b, f);		
			count_XOR_num((*b), cnt, 0, j,R);
		    char str[10];
		    itoa((*b), str, 16);	
		    if(strlen(str) < 2)
		    {
		    	char addzero[10] = "0";
		    	strcat(addzero, str); 
		    	strcpy(str, addzero);
		    }
	    	pb_vec[R][cnt].vec_out += str;
	    	if(j != Matrix_size - 1) pb_vec[R][cnt].vec_out += " ";
	    }
	    pb_vec[R][cnt].vec_out += ") ";
	    pb_vec[R][cnt].XOR_num_out = 4 * (pb_vec[R][cnt].XOR_num_out + 3 * 8);
	    pb_vec[R][cnt].XOR_delay_out = count_XOR_delay();

	    pb_vec[R][cnt].XOR_num_all = pb_vec[R][cnt].XOR_num_in + pb_vec[R][cnt].XOR_num_out;
	    pb_vec[R][cnt].XOR_delay_all = pb_vec[R][cnt].XOR_delay_in + pb_vec[R][cnt].XOR_delay_out ;
	}
	cnt++;
}

int main(){
	SET_BASIC_PARAMETERS;
	precompute_table();

	DWORD a[4] = {0}; a[2] = 1; a[3] = 1; 
	ofstream fout;
	fout.open(txt_name, ios::trunc);
	for(int R=1;R<=255;R++)
		genexpr(R);
	int Max_ele = 256;
	for(int i = 2; i < Max_ele; i++)
	{
		a[0] = i;
		for(int j = i; j < Max_ele; j++)
		{
			a[1] = j;
			inverse_vec(a);
		}
	}
	for(int R=1;R<=255;R++)
	{
		fout<<endl<<endl<<"******************"<<R<<"******************"<<endl<<endl;
		sort(pb_vec[R], pb_vec[R] + cnt, cmp);
		for(int i = 0; i < cnt; i++)
		{
			if(pb_vec[R][i].XOR_num_all < 448 && pb_vec[R][i].XOR_delay_all <= 7)  // 152 + 440
				fout << pb_vec[R][i].vec_in << pb_vec[R][i].XOR_num_in << " " << pb_vec[R][i].XOR_delay_in << " & "
					<< pb_vec[R][i].vec_out << pb_vec[R][i].XOR_num_out << " " << pb_vec[R][i].XOR_delay_out << " & "
					<< pb_vec[R][i].XOR_num_all << " " << pb_vec[R][i].XOR_delay_all << "\n";
		}
	}

	fout.close();
	
	return 0;
}
