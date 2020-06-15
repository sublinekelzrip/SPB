#include<bits/stdc++.h>
using namespace std;
int mul[260][260],exprAB[260][8][8][8],now[260][8][8],ny[260];



void genexpr(int R)
{
	int mod=451;
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
                if(mul[R][mul[1<<i][1<<j]]&(1<<k))
                {
                    exprAB[R][k][j][now[R][k][j]]=i;
                    now[R][k][j]++;
                }
    for(int i=0;i<8;i++)
    {
        for(int j=0;j<8;j++)
        {
            cout<<"{";
            for(int k=0;k<now[R][i][j];k++)
                cout<<exprAB[R][i][j][k]<<' ';
            cout<<"}";
        }
        cout<<endl;
    }
}
vector<int> v[8];
int p;
string s[8];
void count_XOR_num(unsigned long dw,char c,int R)
{
	int arr[8] = {0}; int aaa = dw;
	int w = 0;
	while(dw > 0){
		if(dw & 1) arr[w] = 1;
        else
        arr[w]=0;
		dw >>= 1;
		w++;
	}
	int expr_size = 8, ele_num = 8, XOR_num = 0;
	for(int i = 0; i < expr_size; i++)
	{
        v[i].clear();
		int XOR_num_per_row = 0;
		for(int j = 0; j < expr_size; j++)
		{
			int XOR_num_per_ele = 0;
			for(int k = 0; k < now[R][i][j]; k++)
			{
				if(arr[exprAB[R][i][j][k]] == 1) XOR_num_per_ele ^= 1;
			}
            if(XOR_num_per_ele)
                v[i].push_back(j);
			XOR_num_per_row += XOR_num_per_ele;
		}
        for(int j=0;j<v[i].size();j++)
        {
            s[i]+=c;
            s[i]+="_";
            s[i]+=((char)(v[i][j]+48));
            s[i]+="+";
        }
        cout<<"x^"<<i<<endl;
		XOR_num += XOR_num_per_row;
	}
	XOR_num -= 8;	
    p+=XOR_num;
}
int main()
{
    genexpr(36);
    cout<<ny[111];
    count_XOR_num(0xd6,'a',36);
    count_XOR_num(0x08,'b',36);
    count_XOR_num(0x6f,'c',36);
    count_XOR_num(0xde,'d',36);
    cout<<p<<endl;
    cout<<"\\begin{gather*}\n";
    cout<<"\\begin{split}\n";
    for(int i=7;i>=0;i--)
        cout<<"("<<s[i].substr(0,s[i].length()-1)<<")\\\\"<<endl;
    cout<<"\\end{split}\n";        
    cout<<"\\end{gather*}\n";
}