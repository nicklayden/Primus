#include "gui.hpp"

void DrawParticles(sf::RenderWindow* window, Body_ctr bodies)
{
  // Draws the particles onto the SFML window. No shit.
  // Modelling each particle as a circleshape, and giving them a large
  // enough radius so that you can actually see them on the large window:
  sf::CircleShape particle(6e9); //units are scaled to the size of the viewport.
  for (size_t i = 0; i < bodies.size(); i++) {
    if (bodies[i]->type > 0) {
      particle.setPosition(bodies[i]->rx, bodies[i]->ry);
      if (bodies[i]->name == "basic") {
        particle.setFillColor(sf::Color::Red);
      } else {
        particle.setFillColor(sf::Color::Yellow);
      }
      window->draw(particle);
    }
  }
}
