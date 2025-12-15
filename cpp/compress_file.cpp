// -*- coding: utf-8 -*-
// compile: g++ -O3 compress_file.cpp -o compress_file
// Usage: compress_vels input_file

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <limits>

using namespace std;

// 辅助函数：跳过指定行数
void skip_lines(istream& pStream, size_t pLines) {
    string s;
    for (; pLines; --pLines) {
        if (!getline(pStream, s)) break;
    }
}

int main(int argc, char** argv) {

    // 1. 参数检查
    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " <lammps_dump_file>" << endl;
        return 1;
    }

    string input_filename = argv[1];
    cout << "Reading from file: " << input_filename << endl;

    // 2. 打开输入文件
    ifstream infile(input_filename.c_str());
    if (!infile.is_open()) {
        cerr << "Error: Could not open input file " << input_filename << endl;
        return 1;
    }

    // 3. 获取原子数量 (N)
    // 读取头文件的前几行来解析原子数
    // LAMMPS dump 格式通常如下:
    // ITEM: TIMESTEP
    // ...
    // ITEM: NUMBER OF ATOMS
    // <N>
    // ...
    
    string line;
    int N = 0;
    
    // 跳过前3行 (TIMESTEP行 和 TIMESTEP值 和 "ITEM: NUMBER OF ATOMS")
    skip_lines(infile, 3);
    infile >> N; // 读取第4行：原子数量
    
    if (infile.fail() || N <= 0) {
        cerr << "Error: Failed to read number of atoms or N is invalid." << endl;
        return 1;
    }
    
    cout << "Detected Number of Atoms per frame: " << N << endl;

    // 4. 重置文件指针到开头，准备正式处理
    infile.clear();                 // 清除 EOF 标志
    infile.seekg(0, ios::beg);      // 回到文件头

    // 5. 打开三个二进制输出文件
    ofstream out_vx("vx.bin", ios::binary);
    ofstream out_vy("vy.bin", ios::binary);
    ofstream out_vz("vz.bin", ios::binary);

    if (!out_vx.is_open() || !out_vy.is_open() || !out_vz.is_open()) {
        cerr << "Error: Could not create output binary files." << endl;
        return 1;
    }

    cout << "Writing to vx.bin, vy.bin, vz.bin (binary double precision)..." << endl;

    // 6. 分配缓冲区 (一次性读取一帧的数据再写入，提高IO效率)
    // 使用 double 存储速度
    vector<double> buf_vx(N);
    vector<double> buf_vy(N);
    vector<double> buf_vz(N);

    int id, type;
    double vx, vy, vz;
    size_t frame_count = 0;

    // 7. 主循环：逐帧读取
    while (infile.peek() != EOF) {
        // 尝试读取第一行，如果失败则退出（处理文件末尾的空行）
        // LAMMPS 头有9行
        // 我们这里用 getline 逐行跳过，确保文件流状态正确
        bool header_ok = true;
        for(int i=0; i<9; ++i) {
            if(!getline(infile, line)) {
                header_ok = false; 
                break; 
            }
        }
        if (!header_ok) break; // 文件结束或读取错误

        // 读取 N 个原子的数据
        // 格式: id type vx vy vz
        for (int i = 0; i < N; ++i) {
            infile >> id >> type >> vx >> vy >> vz;
            //infile  >> vx >> vy >> vz;
            
            // 存入缓冲区 (注意：这里假设dump文件已经是 sort id 的)
            // 如果dump文件没有sort，这里存的顺序就是文件里的顺序
            buf_vx[i] = vx;
            buf_vy[i] = vy;
            buf_vz[i] = vz;
        }

        // 读取完一行数据后，流指针停在最后一个数字后，需要跳过这一行的换行符
        // 以便下一轮循环能正确读取 header
        string dummy;
        getline(infile, dummy); 

        // 写入二进制文件
        // reinterpret_cast 将 double* 转换为 char* 供 write 使用
        // sizeof(double) * N 表示写入的总字节数
        out_vx.write(reinterpret_cast<char*>(buf_vx.data()), N * sizeof(double));
        out_vy.write(reinterpret_cast<char*>(buf_vy.data()), N * sizeof(double));
        out_vz.write(reinterpret_cast<char*>(buf_vz.data()), N * sizeof(double));

        frame_count++;
        if (frame_count % 100 == 0) {
            cout << "Processed frame: " << frame_count << "\r" << flush;
        }
    }

    cout << "\nDone! Total frames processed: " << frame_count << endl;

    // 8. 关闭文件
    infile.close();
    out_vx.close();
    out_vy.close();
    out_vz.close();

    return 0;
}
