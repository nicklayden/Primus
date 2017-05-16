#ifndef GUI_HPP
#define GUI_HPP
#include <SFML/Graphics.hpp>
#include <vector>
#include <sstream>
#include "Particle.hpp"
#include "simulator.hpp"

typedef std::vector<Particle*> Body_ctr;

// Functions for drawing particles and creating the SFML window
void DrawParticles(sf::RenderWindow* window, Body_ctr& bodies);
void DrawBoxes(sf::RenderWindow* window, std::vector<std::vector<double> > nodes);


template<class T>
std::string NumberToString(std::string text, T Number, std::string unit)
{
  std::ostringstream ss;
  ss << text << " " << Number << " " << unit;
  return ss.str();
}

class SFMLDisplay
{
public:
  SFMLDisplay(int wWidth, int wHeight, Simulator* sim);
  ~SFMLDisplay();

  void draw(Body_ctr& bodies);


private:
  Simulator* sim;
  sf::RenderWindow mainWindow;
  sf::RenderWindow statsWindow;
  double particleSize = 6e9;
};


#endif /* end of include guard: GUI_HPP */
