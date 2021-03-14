
//Compile with : g++  *.cpp -IPathToSDL2\include -LPathToSDL2\lib -Wall -lmingw32 -lSDL2main -lSDL2 -o main
//Example: g++  *.cpp -IC:\SDL2\include -LC:\SDL2\lib -Wall -lmingw32 -lSDL2main -lSDL2 -o main
// #include <iostream>
#include <bits/stdc++.h>
#include <SDL2/SDL.h>
#include "Fluid.h"
#define endl "\n"
#define WIDTH 640
#define HEIGHT 640
int main(int argc, char *argv[]){
  srand(time(nullptr));
  if(SDL_Init(SDL_INIT_EVERYTHING)){
    std::cout << "Error : " << SDL_GetError() << endl;
    return 0;
  }
  SDL_Window *wind = SDL_CreateWindow("Fluid Grid", SDL_WINDOWPOS_CENTERED,SDL_WINDOWPOS_CENTERED,WIDTH,HEIGHT,0);
  if (!wind)
  {
    printf("Error creating window: %s\n", SDL_GetError());
    SDL_Quit();
    return 0;
  }
  uint32_t flags = SDL_RENDERER_ACCELERATED;
  SDL_Renderer *render = SDL_CreateRenderer(wind,-1,flags);
  if (!render)
  {
    printf("Error creating renderer: %s\n", SDL_GetError());
    SDL_DestroyWindow(wind);
    SDL_Quit();
    return 0;
  }
  bool running = true, pressed = false;
  SDL_Event event;
  int N = 64;              // Size of the fluid cell N*N
  float diffusion = 0;   // Diffusion constant
  float viscosity = 0.000001;     // Viscosity constant
  float dt = 0.1;     // Timestep
  int Scale = 10;
  Fluid fluid = Fluid(N,diffusion,viscosity,dt);
  SDL_Rect rect = {(WIDTH-Scale)/2 , (HEIGHT-Scale)/2, Scale,Scale}; //Square object responsible for displaying the density of each cells
  float amtX, amtY;
  int mouseX, mouseY;
  while(running){
    while(SDL_PollEvent(&event)){
      switch(event.type){
        case SDL_MOUSEMOTION:
          SDL_GetMouseState(&mouseX,&mouseY);
          // If you want to add the fluid by clicking at your cursor position instead of following the cursor,
          // Comment these 3 lines and uncomment lines 70 - 74 
          amtX = amtY = rand()%11-5;
          fluid.addDensity(mouseX/Scale,mouseY/Scale, 100);
          fluid.addVelocity(mouseX/Scale,mouseY/Scale,amtX,amtY);
          break;
        case SDL_MOUSEBUTTONDOWN:
          if(event.button.button == SDL_BUTTON_LEFT){
            pressed = true;
          }
          break;
        case SDL_MOUSEBUTTONUP:
          pressed = false;
          break;
        case SDL_QUIT:
          running = false;
          break;
      }
    }
    SDL_SetRenderDrawColor(render,0,0,0,255);
    SDL_RenderClear(render);
    // if(pressed){
    //   amtX = amtY = rand()%11-5;
    //   fluid.addDensity(mouseX/Scale,mouseY/Scale, 100);
    //   fluid.addVelocity(mouseX/Scale,mouseY/Scale,amtX,amtY);
    // }
    fluid.step();
    for(int i = 0; i < fluid.Size(); ++i){
      for(int j = 0; j< fluid.Size(); ++j){
        int x = j*Scale;
        int y = i*Scale;
        float density = fluid.Density(j,i);
        rect.x = x;
        rect.y = y;
        SDL_SetRenderDrawBlendMode(render,SDL_BLENDMODE_BLEND);
        SDL_SetRenderDrawColor(render, 100, 255, 255, density);
        SDL_RenderFillRect(render,&rect);
      }
    }
    /* Draw to window and loop */
    SDL_RenderPresent(render);
  }
  
  SDL_Quit();
  return 0;
}
