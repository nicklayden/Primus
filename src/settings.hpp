#ifndef SETTINGS_HPP
#define SETTINGS_HPP
#pragma once

#include <string>
#include <fstream>
#include <iostream>

using std::string;

template<class T>
class setting
{
    public:
        setting(string name, T val)
        :name(name), val(val)
        {}

        string name;
        T val;

        friend std::ostream& operator<<(std::ostream& out, const setting& s)
        {
            return out << s.name << " " << s.val << std::endl;
        }

        setting& operator=(T rhs)
        {
            val = rhs;
        }

        setting& operator*(setting& rhs)
        {
            return val*rhs.val;
        }

};




#endif /* end of include guard: SETTINGS_HPP */
