


#include "HX711.h"


#define potpin A2

int set_speed,motor_speed,disp_speed;

HX711 scale;
float pos=0;

unsigned long t ;
unsigned long tGo ;

float stepCount = 0;

int steppersteps= 200;

const int buttonGo = 8;//buttons to move the motor, now analog
long previousMillis = 0;
long currentMillis = 0;


//Calibrating the load cells
float ReadingA_Strain1 = 36600;
float LoadA_Strain1 = 20; // (Kg,lbs..)
float ReadingB_Strain1 = 185000;
float LoadB_Strain1 = 100; // (Kg,lbs..)

float ReadingA_Strain2 = 37;
float LoadA_Strain2 = 0.5094; // (Kg,lbs..)
float ReadingB_Strain2 = 54;
float LoadB_Strain2 = 1.0094; // (Kg,lbs..)

void setup() {
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

   scale.begin(2, 3);
  scale.tare();
  


 
  pinMode(buttonGo, INPUT);
  
  digitalWrite(buttonGo, HIGH);//enable internal pullups



  scale.tare();



}



void loop() {



  boolean go = false;
  
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
     tGo = millis();
      
  }
  else {
    pressed = false;

  }




  while (pressed == true  ) {
   
    unsigned long currentMillis = millis();
   
    t = currentMillis-tGo;

     

        pos =0; //pos in mm

         Serial.print(t);
        Serial.print(",");
         Serial.print(pos,3);
        Serial.print(",");
        Serial.print(steppersteps);
        Serial.print(",");
        
        
        float newReading_Strain1 = scale.get_units(1);
  float load_Strain1 = ((LoadB_Strain1 - LoadA_Strain1) / (ReadingB_Strain1 - ReadingA_Strain1)) *
                       (newReading_Strain1 - ReadingA_Strain1) + LoadA_Strain1;
   Serial.print( newReading_Strain1);
        
Serial.print(",");
    Serial.println( load_Strain1,3);
        
 
       
        stepCount=0;
        



   }
     

}
