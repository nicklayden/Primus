/*
  Display window for Primus.
*/
#include "display.hpp"
#include "simulator.hpp"


Display::Display(Simulator* simulation, sf::Font font)
:simulation(simulation), font(font), window(sf::VideoMode(800,800), "test")
{

  std::cout << "Display window created." << std::endl;
  if (!font.loadFromFile("OpenSans-Regular.ttf")) {
    std::cout << "Failed getting text font." << std::endl;
  }
}

Display::Display(sf::Font font)
:font(font), window(sf::VideoMode(800,800), "test")
{
  std::cout << "Display window created." << std::endl;
  std::cout << "Basic constructor called." << std::endl;
}

Display::~Display()
{
  std::cout << "Window destroyed." << std::endl;
}

void Display::MainLoop()
{
  sf::Event event;
  ProcessEvents(event);
  window.clear(backgroundColour);
  // Draw commands here


  DrawText();
  // end draw commands
  window.display();
}

void Display::ProcessEvents(sf::Event event)
{
  while (window.pollEvent(event)){
    if (event.type == sf::Event::Closed)
    {
      window.close();
    } else if (event.type == sf::Event::KeyPressed) {
      // Handle keyboard presses
        // A key
      if (event.key.code == sf::Keyboard::A && drawNodes == true) {
        drawNodes = false;
      } else if (event.key.code == sf::Keyboard::A && drawNodes == false) {
        drawNodes = true;
        // S key
      } else if (event.key.code == sf::Keyboard::S && showStats == true) {
        showStats = false;
      } else if (event.key.code == sf::Keyboard::S && showStats == false) {
        showStats = true;
      }
    }
  }
}

void Display::DrawText()
{
  // Draw the statistics from the simulation on the screen.
  timer.setFont(font);
  timer.setString("0");
  timer.setCharacterSize(12);
  timer.setColor(sf::Color::White);
  timer.setPosition(0,12);

  simspeed.setFont(font);
  simspeed.setString("0");
  simspeed.setCharacterSize(12);
  simspeed.setColor(sf::Color::White);
  simspeed.setPosition(0,24);

  numBodies.setFont(font);
  numBodies.setString("0");
  numBodies.setCharacterSize(12);
  numBodies.setColor(sf::Color::White);
  numBodies.setPosition(0,36);

  forcesperstep.setFont(font);
  forcesperstep.setString("0");
  forcesperstep.setCharacterSize(12);
  forcesperstep.setColor(sf::Color::White);
  forcesperstep.setPosition(0,48);

  flopsTimer.setFont(font);
  flopsTimer.setString("0");
  flopsTimer.setCharacterSize(12);
  flopsTimer.setColor(sf::Color::White);
  flopsTimer.setPosition(0,60);

  currentTimestep.setFont(font);
  currentTimestep.setString("0");
  currentTimestep.setCharacterSize(12);
  currentTimestep.setColor(sf::Color::White);
  currentTimestep.setPosition(0,72);

  timeChanges.setFont(font);
  timeChanges.setString("0");
  timeChanges.setCharacterSize(12);
  timeChanges.setColor(sf::Color::White);

  window.draw(currentTimestep);
  window.draw(timeChanges);
  window.draw(timer);
  window.draw(numBodies);
  window.draw(simspeed);
  window.draw(forcesperstep);
  window.draw(flopsTimer);

}