#ifndef GUI_HPP
#define GUI_HPP
#include <SFML/Graphics.hpp>
#include <vector>
#include <sstream>
#include <iostream>
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


/******************************************************************************/
//    class:  SFMLDisplay
//    purpose: Draws and displays particle simulation in real time.
//
class SFMLDisplay
{
public:
  SFMLDisplay(Simulator* sim, sf::Font* font);
  ~SFMLDisplay();

  // Public interface for drawing the window.
  void drawParticles();
  void drawText();
  void displayLoop();
  void displaySingleFrame();


private:
  // Private member functions
  void handleEvents(sf::Event event);
  void initStats();

  // Private member variables
  float particleSize = 6e9;
  Simulator* sim;
  sf::RenderWindow mainWindow;
  sf::View mainView;
  sf::View statView;
  sf::Font* font;
  sf::Color textColor = sf::Color::White;
  std::vector<sf::Text*> simStats;
  int fontSize = 12;
  // Keyboard switch flags
  bool displayStats = true;
  bool displayNodes = true;

};


#endif /* end of include guard: GUI_HPP */
