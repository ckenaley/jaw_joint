#include <AccelStepper.h>

//1201

unsigned long JL = 30.09;

int adjust = 500;

int maxSpeed = 300; 
//300w
//50

int gosteps= (tan(0.349)*JL)*(1000/6.3); //1000 = 6.3mm

int count = 1   ;//number of cycles
int i = 1;
const int buttonGo = 11;//buttons to move the motor, now analog
const int buttonPush = 12;//buttons to move the motor
const int buttonPull = 13;//buttons to move the motor

unsigned long t ;
unsigned long tGo ;


long previousMillis = 0;
long currentMillis = 0;

// Define a stepper and the pins it will use
AccelStepper stepper(AccelStepper::DRIVER, 8, 7);

int pos = 300;

//500 steps is 3 mm
//Calibrating the load cells
float ReadingA_Strain1 = -98900;
float LoadA_Strain1 = 50; // (Kg,lbs..)
float ReadingB_Strain1 = -197800;
float LoadB_Strain1 = 100; // (Kg,lbs..)


void setup()
{  


Serial.begin(115200);
      Serial.print("ms");
        Serial.print(",");
      Serial.print("pos");
        Serial.print(",");
        Serial.print("steps");
        Serial.print(",");
        
      
    Serial.print( "strain" );
     Serial.print(",");
      Serial.println( "g" );

  
 
  pinMode(buttonGo, INPUT);
  pinMode(buttonPush, INPUT);
  pinMode(buttonPull, INPUT);
  digitalWrite(buttonGo, HIGH);//enable internal pullups
  digitalWrite(buttonPull, HIGH);//enable internal pullups
  digitalWrite(buttonPush, HIGH);//enable internal pullups
}

void loop()
{

  boolean go = false;
  boolean push = false;
  boolean pull = false;

  /////// go  button signal
  if (digitalRead(buttonGo) == LOW) {
    go = true;
   
  }
  else {
    go = false;
  }

  boolean pressed = false;
  if (go == true) {
    pressed = true;
     tGo=millis();
  }
  else {
    pressed = false;
    
  }

   /////// pull  button signal
  if (digitalRead(buttonPull) == LOW) {
    pull = true;
  }
  else {
    pull = false;
  }

  while (pull == true) {
    stepper.setMaxSpeed(1000);
    stepper.setMaxSpeed(1000);
  stepper.setAcceleration(5000);
 stepper.setCurrentPosition(0);
   stepper.moveTo(adjust);
   stepper.run();

   
      
    
   // delay(300);

    pull = false;
   
  }

   /////// push  button signal
  if (digitalRead(buttonPush) == LOW) {
    push = true;
    
  }
  else {
    push = false;
  }

  while (push == true) {
    stepper.setMaxSpeed(1000);
  stepper.setAcceleration(5000);
 stepper.setCurrentPosition(0);
   stepper.moveTo(-adjust);
   stepper.run();
    
      
    
   // delay(300);

    push = false;
   
  }

  while (pressed == true) {
     if (i <= count ) {
    stepper.setMaxSpeed(maxSpeed);
  stepper.setAcceleration(20000);
      stepper.setCurrentPosition(0);
        
stepper.runToNewPosition(gosteps);

    stepper.stop();

    unsigned long currentMillis = millis();
   
    t = currentMillis-tGo;
    //Serial.println(t);
    delay(100);
  stepper.runToNewPosition(0); // Cause an overshoot then back to 0
  
  }
  
      i++;
      stepper.stop();
      delay(1001);
     if(i>count){
       pressed=false;
     }
     }


}
