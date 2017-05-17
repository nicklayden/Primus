#include "gui.hpp"
#include "AstroConstants.hpp"
using namespace constants;

void DrawParticles(sf::RenderWindow* window, Body_ctr& bodies)
{
  // Draws the particles onto the SFML window. No shit.
  // Modelling each particle as a circleshape, and giving them a large
  // enough radius so that you can actually see them on the large window:
  sf::CircleShape particle(0.04); //units are scaled to the size of the viewport.
  for (size_t i = 0; i < bodies.size(); i++) {
    if (bodies[i]->type > 0) {
      particle.setPosition(bodies[i]->rx/au, bodies[i]->ry/au);
      if (bodies[i]->name == "basic") {
        particle.setFillColor(sf::Color::Red);
      } else {
        particle.setFillColor(sf::Color::Yellow);
      }
      window->draw(particle);
    }
  }
}

void DrawBoxes(sf::RenderWindow* window, std::vector<std::vector<double> > nodes)
{
  for (size_t i = 0; i < nodes.size(); i++) {
    sf::RectangleShape rect(sf::Vector2f(nodes[i][2]/au,nodes[i][3]/au));
    rect.setPosition(nodes[i][0]/au,nodes[i][1]/au);
    rect.setOutlineColor(sf::Color::Green);
    rect.setOutlineThickness(-0.02);
    rect.setFillColor(sf::Color::Transparent);
    window->draw(rect);
  }
}


SFMLDisplay::SFMLDisplay(int wWidth, int wHeight, Simulator* sim)
:mainWindow(sf::VideoMode(wWidth, wHeight), "Simulation"), sim(sim),
 statsWindow(sf::VideoMode(500,150),"Simulation Stats")
{
}

void SFMLDisplay::draw(Body_ctr& bodies)
{
  sf::CircleShape particle(particleSize); //units are scaled to the size of the viewport.
  for (size_t i = 0; i < bodies.size(); i++) {
    if (bodies[i]->type > 0) {
      particle.setPosition(bodies[i]->rx, bodies[i]->ry);
      if (bodies[i]->name == "basic") {
        particle.setFillColor(sf::Color::Red);
      } else {
        particle.setFillColor(sf::Color::Yellow);
      }
      mainWindow.draw(particle);
    }
  }
}
