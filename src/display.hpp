/*
    Display window for Primus. This class handles all SFML display and drawing.
*/

#include <SFML/Graphics.hpp>
#include "simulator.hpp"

class Display
{
  public:
    Display(Simulator* simulation, sf::Font font);
    Display(sf::Font font, double w, double h);
    ~Display();
    void MainLoop(std::vector<Particle*> bodies);

    double particleSize = 6e9;

    bool showText = false;
    bool showStats = false;
    bool drawNodes = false;

    sf::RenderWindow window;
    sf::Text title;
    sf::Color backgroundColour = sf::Color::Black;
    sf::View mainview;
    sf::View statview;

  private:
    Simulator* simulation;
    void DrawSimulation();
    void DrawText();
    void Scatter();
    void GetBoxDimensions();
    void EnterNodeTree();
    void ProcessEvents(sf::Event event);
    void DrawParticles(std::vector<Particle*>& bodies);


    int charsize = 12;

    sf::Font font;
    sf::Text timer;
    sf::Text simspeed;
    sf::Text numBodies;
    sf::Text forcesperstep;
    sf::Text flopsTimer;
    sf::Text currentTimestep;
    sf::Text timeChanges;

};