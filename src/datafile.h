#ifndef _DATAFILE_H
#define _DATAFILE_H

#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <map>

struct datafile
{
    int n_skipped_line;
    int n_column;
    
    std::vector <std::string> skipped_line;
    std::vector <std::vector<std::string> > data;

    void load(char *fn, int skip=1)
    {
        n_skipped_line=skip;
        std::ifstream in(fn);
        char buf[0x1000];
        std::string tempstr;

        // read skipped line
        for(int i=0;i<n_skipped_line;i++)
        {
            in.getline(buf,0x1000);
            tempstr=buf;
            skipped_line.push_back(tempstr);
        }

        // test read first raw and get column count
        std::streampos spos=in.tellg();
        in.getline(buf,0x1000);
        int slen=strlen(buf);
        n_column=1;
        for(int i=0;i<slen-1;i++) if(buf[i]<=32&&buf[i+1]>32) n_column++;
        in.seekg(spos);

        // initialize data
        data.resize(n_column);
        
        // load data
        while(!in.eof())
        {
            for(int i=0;i<n_column;i++)
            {
                in>>tempstr;
                if(in.eof()) break;
                data[i].push_back(tempstr);
            }
        }
        in.close();
    };

    void get_data(std::vector<float> &v, int column)
    {
        int size=data[column].size();
        v.resize(size);
        for(int i=0;i<size;i++) v[i]=atof(data[column][i].c_str());
    };

    void get_data(std::vector<int> &v, int column)
    {
        int size=data[column].size();
        v.resize(size);
        for(int i=0;i<size;i++) v[i]=atoi(data[column][i].c_str());
    };
    
    void get_data(std::vector<std::string> &v, int column)
    {
        int size=data[column].size();
        v.resize(size);
        for(int i=0;i<size;i++) v[i]=data[column][i];
    };

    template <class T>
    void set_map(std::map<std::string,T> &m, int key_column, int data_column)
    {
        std::vector<T> d;
        get_data(d,data_column);
        int size=data[key_column].size();
        for(int i=0;i<size;i++) m[data[key_column][i]]=d[i];
    };
};
#endif
