#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <functional>
using namespace std;

// 質点の情報を保持するためのクラス
class MyPoint {
    public:
        double m, x, y, z, vx, vy, vz;
        MyPoint(double m, double x, double y, double z, double vx, double vy, double vz){
            this->m = m;
            this->x = x;
            this->y = y;
            this->z = z;
            this->vx = vx;
            this->vy = vy;
            this->vz = vz;
        }
};

// 連立方程式のパラメータを保持するクラス
class Var {
    public:
        vector<double> values;
        double init;
        function<double(vector<double>)> func;
        Var( double init, function<double(vector<double>)> f ) {
            this->values.push_back(init);
            this->init = init;
            this->func = f;
        }
        void clear() {
            this->values.clear();
            this->values.push_back( this->init );
        }
};

vector<string> str_split( string target, string sep = " ");
void runge_kutta_2( vector<Var>& vars, double h, double t_end, vector<MyPoint>& points );
void runge_kutta_4( vector<Var>& vars, double h, double t_end, vector<MyPoint>& points );
double get_energy( vector<MyPoint>& points, vector<Var>& vars, int i );

int main(){
    // ifstream ifs("./report3_1.dat");
    ifstream ifs("./report3_2.dat");
    string line1;
    vector<string> info, lines; // Variables for reading dat-file.
    vector<MyPoint> points;


    if ( ifs.fail() ) return 0;

    getline( ifs, line1 );
    info = str_split(line1);

    int num = stoi(info.at(0));
    double dt = stod(info.at(1));
    double t_end = stod(info.at(2));

    for ( int i = 0; i < num; i++ ) {
        string tmp_s;
        getline( ifs, tmp_s );
        vector<string> tmp_d = str_split(tmp_s);
        points.push_back( *(new MyPoint(stod(tmp_d.at(0)), stod(tmp_d.at(1)), stod(tmp_d.at(2)), stod(tmp_d.at(3)), stod(tmp_d.at(4)), stod(tmp_d.at(5)), stod(tmp_d.at(6)))) );
    }

    vector<Var> vars;
    const int size = points.size();
    for ( int i = 0; i < size; i++ ) {
        vars.push_back( *(new Var(points.at(i).x, [i](vector<double> params){
            return params.at(6*i+3);
        })) );
        vars.push_back( *(new Var(points.at(i).y, [i](vector<double> params){
            return params.at(6*i+4);
        })) );
        vars.push_back( *(new Var(points.at(i).z, [i](vector<double> params){
            return params.at(6*i+5);
        })) );
        vars.push_back( *(new Var(points.at(i).vx, [i, points, size](vector<double> params){
            double tmp = 0;
            for ( int j = 0; j < size; j++ ) {
                if ( i == j ) continue;
                double dx = params.at(6*j)-params.at(6*i);
                double dy = params.at(6*j+1)-params.at(6*i+1);
                double dz = params.at(6*j+2)-params.at(6*i+2);
                tmp += points.at(j).m * dx / pow(dx*dx + dy*dy + dz*dz, 1.5);
            }
            return tmp;
        })) );
        vars.push_back( *(new Var(points.at(i).vy, [i, points, size](vector<double> params){
            double tmp = 0;
            for ( int j = 0; j < size; j++ ) {
                if ( i == j ) continue;
                double dx = params.at(6*j)-params.at(6*i);
                double dy = params.at(6*j+1)-params.at(6*i+1);
                double dz = params.at(6*j+2)-params.at(6*i+2);
                tmp += points.at(j).m * dy / pow(dx*dx + dy*dy + dz*dz, 1.5);
            }
            return tmp;
        })) );
        vars.push_back( *(new Var(points.at(i).vz, [i, points, size](vector<double> params){
            double tmp = 0;
            for ( int j = 0; j < size; j++ ) {
                if ( i == j ) continue;
                double dx = params.at(6*j)-params.at(6*i);
                double dy = params.at(6*j+1)-params.at(6*i+1);
                double dz = params.at(6*j+2)-params.at(6*i+2);
                tmp += points.at(j).m * dz / pow(dx*dx + dy*dy + dz*dz, 1.5);
            }
            return tmp;
        })) );
    }

    // 2nd Runge-Kutta
    runge_kutta_2( vars, dt, t_end, points );

    // 4th Runge-Kutta
    runge_kutta_4( vars, dt, t_end, points );
    return 0;
}

vector<string> str_split( string target, string sep ) {
    vector<string> list;
    size_t offset = string::size_type(0);
    while ( true ) {
        int pos = target.find(sep, offset);
        if (pos == string::npos) {
            list.push_back(target.substr(offset));
            break;
        }
        if ( target.substr(offset, pos - offset) != "" ) {
            list.push_back(target.substr(offset, pos - offset));
        }
        offset = pos + sep.length();
    }
    return list;
}

void runge_kutta_2( vector<Var>& vars, double h, double t_end, vector<MyPoint>& points ) {
    const int size = vars.size();
    const double init_energy = get_energy( points, vars, 0 );

    for ( double t = h; t <= t_end; t += h ) {
        vector<double> now_vals;
        for ( int i = 0; i<size; i++ ) {
            now_vals.push_back( vars.at(i).values.back() );
        }

        vector<double> k1;
        for ( int i = 0; i < size; i++ ) {
            k1.push_back( h * vars.at(i).func(now_vals) );
        }

        vector<double> tmp;
        for ( int i = 0; i < size; i++ ) {
            tmp.push_back( now_vals.at(i) + k1.at(i) );
        }

        vector<double> k2;
        for ( int i = 0; i < size; i++ ) {
            k2.push_back( h * vars.at(i).func(tmp) );
        }

        for ( int i = 0; i < size; i++ ) {
            vars.at(i).values.push_back( vars.at(i).values.back() + ( k1.at(i) + k2.at(i) ) / 2.0 );
        }
    }

    ofstream ofs("./runge_kutta_2.dat");
    ofs << "# 2次のルンゲ・クッタ法" << endl;
    for ( int i = 0; i < t_end/h; i++ ) {
        ofs << i*h << flush;
        for ( int j = 0; j < size/6; j++ ) {
            ofs << " " << vars.at(6*j).values.at(i) << flush;
            ofs << " " << vars.at(6*j+1).values.at(i) << flush;
            ofs << " " << vars.at(6*j+2).values.at(i) << flush;
        }
        ofs << " " << get_energy(points, vars, i) - init_energy << endl;
    }

    for ( int i = 0; i < size; i++ ) {
        vars.at(i).clear();
    }
}

void runge_kutta_4( vector<Var>& vars, double h, double t_end, vector<MyPoint>& points ) {
    const int size = vars.size();
    const double init_energy = get_energy( points, vars, 0 );

    for ( double t = h; t <= t_end; t += h ) {
        vector<double> now_vals;
        for ( int i = 0; i<size; i++ ) {
            now_vals.push_back( vars.at(i).values.back() );
        }

        vector<double> k1;
        for ( int i = 0; i < size; i++ ) {
            k1.push_back( h * vars.at(i).func(now_vals) );
        }

        vector<double> tmp;
        for ( int i = 0; i < size; i++ ) {
            tmp.push_back( now_vals.at(i) + k1.at(i)/2 );
        }

        vector<double> k2;
        for ( int i = 0; i < size; i++ ) {
            k2.push_back( h * vars.at(i).func(tmp) );
        }

        tmp.clear();
        for ( int i = 0; i < size; i++ ) {
            tmp.push_back( now_vals.at(i) + k2.at(i)/2 );
        }

        vector<double> k3;
        for ( int i = 0; i < size; i++ ) {
            k3.push_back( h * vars.at(i).func(tmp) );
        }

        tmp.clear();
        for ( int i = 0; i < size; i++ ) {
            tmp.push_back( now_vals.at(i) + k3.at(i) );
        }

        vector<double> k4;
        for ( int i = 0; i < size; i++ ) {
            k4.push_back( h * vars.at(i).func(tmp) );
        }

        for ( int i = 0; i < size; i++ ) {
            vars.at(i).values.push_back( vars.at(i).values.back() + ( k1.at(i) + k2.at(i)*2 + k3.at(i)*2 + k4.at(i) ) / 6.0 );
        }
    }

    double en_error;
    ofstream ofs("./runge_kutta_4.dat");
    ofs << "# 4次のルンゲ・クッタ法" << endl;
    for ( int i = 0; i < t_end/h; i++ ) {
        ofs << i*h << flush;
        en_error = 0.0;
        for ( int j = 0; j < size/6; j++ ) {
            ofs << " " << vars.at(6*j).values.at(i) << flush;
            ofs << " " << vars.at(6*j+1).values.at(i) << flush;
            ofs << " " << vars.at(6*j+2).values.at(i) << flush;
            en_error += 0.5 * points.at(j).m * (pow(vars.at(6*j+3).values.at(i), 2)+pow(vars.at(6*j+4).values.at(i), 2)+pow(vars.at(6*j+5).values.at(i), 2));
        }
        ofs << " " << get_energy(points, vars, i) - init_energy << endl;
    }

    for ( int i = 0; i < size; i++ ) {
        vars.at(i).clear();
    }
}

double get_energy( vector<MyPoint>& points, vector<Var>& vars, int i ) {
    double ret = 0.0;
    const int size = points.size();

    for ( int j = 0; j < size; j++ ) {
        ret += points.at(j).m * ( pow(vars.at(6*j+3).values.at(i), 2) + pow(vars.at(6*j+4).values.at(i), 2) + pow(vars.at(6*j+5).values.at(i), 2) / 2 );
        for ( int k = 0; k < j; k++ ) {
            double r = pow( pow(vars.at(6*j).values.at(i)-vars.at(6*k).values.at(i), 2) + pow(vars.at(6*j+1).values.at(i)-vars.at(6*k+1).values.at(i), 2) + pow(vars.at(6*j+2).values.at(i)-vars.at(6*k+2).values.at(i), 2), 0.5 );
            ret += points.at(j).m * points.at(k).m / r;
        }
    }
    return ret;
}
