#include "gui.hpp"
#include <string>
#include "AstroConstants.hpp"
using namespace constants;

void DrawParticles(sf::RenderWindow* window, Body_ctr& bodies)
{
  // Draws the particles onto the SFML window. No shit.
  // Modelling each particle as a circleshape, and giving them a large
  // enough radius so that you can actually see them on the large window:
  sf::CircleShape particle(6e9); //units are scaled to the size of the viewport.
  for (size_t i = 0; i < bodies.size(); i++) {
    if (bodies[i]->type > 0) {
      particle.setPosition(bodies[i]->rx, bodies[i]->ry);
      if (bodies[i]->name == "basic") {
        particle.setFillColor(sf::Color::White);
      } else {
        particle.setFillColor(sf::Color::Yellow);
      }
      window->draw(particle);
    }
  }
}

void DrawBoxes(sf::RenderWindow* window, std::vector<std::vector<double> > nodes)
{
  // Draw the boxes in the BH tree.
  double xmax, xmin, ymax, ymin;
  for (size_t i = 0; i < nodes.size(); i++) {
    xmax = (nodes[i][0] + nodes[i][2]);
    xmin = nodes[i][0];
    ymax = (nodes[i][1] + nodes[i][3]);
    ymin = nodes[i][1];
    sf::VertexArray boundingbox(sf::LinesStrip, 5);
    boundingbox[0].position = sf::Vector2f(xmin, ymax);
    boundingbox[1].position = sf::Vector2f(xmax, ymax);
    boundingbox[2].position = sf::Vector2f(xmax, ymin);
    boundingbox[3].position = sf::Vector2f(xmin, ymin);
    boundingbox[4].position = sf::Vector2f(xmin, ymax);

    for (size_t i = 0; i < 5; i++) {
      boundingbox[i].color = sf::Color::Red;
    }
    window->draw(boundingbox);
  }
}

//
// SFMLDisplay::SFMLDisplay(int wWidth, int wHeight, Simulator* sim)
// :mainWindow(sf::VideoMode(wWidth, wHeight), "Simulation"), sim(sim),
//  statsWindow(sf::VideoMode(500,150),"Simulation Stats")
// {
// }

// void SFMLDisplay::draw(Body_ctr& bodies)
// {
//   sf::CircleShape particle(particleSize); //units are scaled to the size of the viewport.
//   for (size_t i = 0; i < bodies.size(); i++) {
//     if (bodies[i]->type > 0) {
//       particle.setPosition(bodies[i]->rx, bodies[i]->ry);
//       if (bodies[i]->name == "basic") {
//         particle.setFillColor(sf::Color::Red);
//       } else {
//         particle.setFillColor(sf::Color::Yellow);
//       }
//       mainWindow.draw(particle);
//     }
//   }
// }
/******************************************************************************/
//    class:  SFMLDisplay
//    purpose: Draws and displays particle simulation in real time.
//

SFMLDisplay::SFMLDisplay(Simulator* sim, sf::Font* font)
:mainWindow(sf::VideoMode(1080,1080),"Simulation"), sim(sim), font(font)
{
  // statView.setCenter(0,0);
  // statView.setSize(100,-100);
  mainView.setCenter(0,0);
  mainView.setSize(sim->w, sim->h);
  mainWindow.setView(mainView);
  initStats();
  displayLoop();
}

void SFMLDisplay::drawParticles()
{
  // Draws the particles in simulation.
  sf::CircleShape particle(particleSize); //units are scaled to the size of the viewport.
  for (size_t i = 0; i < sim->bodies.size(); i++) {
    if (sim->bodies[i]->type > 0) {
      particle.setPosition(sim->bodies[i]->rx, sim->bodies[i]->ry);
      if (sim->bodies[i]->name == "basic") {
        particle.setFillColor(sf::Color::Red);
      } else {
        particle.setFillColor(sf::Color::Yellow);
      }
      mainWindow.draw(particle);
    }
  }
}

void SFMLDisplay::handleEvents(sf::Event event)
{
  while (mainWindow.pollEvent(event)) {
    if (event.type == sf::Event::Closed) {
      mainWindow.close();
    }
    else if (event.type == sf::Event::KeyPressed) {
      if (event.key.code == sf::Keyboard::A && displayNodes == true) {
        displayNodes = false;
      }
      else if (event.key.code == sf::Keyboard::A && displayNodes == false) {
        displayNodes = true;
      }
      else if (event.key.code == sf::Keyboard::S && displayStats == true) {
        displayStats = false;
      }
      else if (event.key.code == sf::Keyboard::S && displayStats == false) {
        displayStats = true;
      }
    }
  } // end while poll
}

void SFMLDisplay::displayLoop()
{
  while (mainWindow.isOpen()) {
    sf::Event event;
    handleEvents(event);
    mainWindow.clear(sf::Color::Black);
    mainWindow.setView(mainView);
    drawParticles();
    if (displayStats) {
      mainWindow.setView(statView);
      drawText();
    }

    mainWindow.display();
  }
}

void SFMLDisplay::drawText()
{
  for (size_t i = 0; i < simStats.size(); i++) {
    mainWindow.draw(*simStats[i]);
  }
}

void SFMLDisplay::initStats()
{
  // Initialize the text displayed for statistics.
  sf::Text* numParticles = new sf::Text;
  numParticles->setFont(*font);
  numParticles->setString(std::to_string(sim->bodies.size()));
  numParticles->setCharacterSize(fontSize);
  numParticles->setPosition(0,0);
  numParticles->setColor(textColor);
  simStats.push_back(numParticles);
}

SFMLDisplay::~SFMLDisplay()
{
  for (size_t i = 0; i < simStats.size(); i++) {
    delete simStats[i];
  }
  std::cout << "Deleted display window." << '\n';
}
