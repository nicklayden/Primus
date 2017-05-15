#ifndef GUI_HPP
#define GUI_HPP
#include <SFML/Graphics.hpp>
#include <vector>
#include <sstream>
#include "Particle.hpp"


typedef std::vector<Particle*> Body_ctr;

// Functions for drawing particles and creating the SFML window
void DrawParticles(sf::RenderWindow* window, Body_ctr bodies);

template<class T>
std::string NumberToString(std::string text, T Number, std::string unit)
{
  std::ostringstream ss;
  ss << text << " " << Number << " " << unit;
  return ss.str();
}



#endif /* end of include guard: GUI_HPP */
