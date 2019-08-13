// declare variables
const int clk = 1;//pin for clock signal
const int data = 2;//pin that sends settings to attenuators
const int latch = 3;//pin for latch. Change this number to change which attenuator is set
int count = 0;
int attenSet[6] = {1,0,0,0,0,0};//settings for attenuation. See zx76-31a-sns data sheet (in datasheet folder). Formatted as [+/- 16dB, +/- 8dB, +/- 4dB, +/- 2dB, +/- 1dB, +/- 0.5dB]
int del = 50;

void setup() {

  //setup pins
  pinMode(clk,OUTPUT);
  pinMode(data,OUTPUT);
  pinMode(LE1,OUTPUT);
  pinMode(LE2,OUTPUT);
  pinMode(LE3,OUTPUT);
  pinMode(LE4,OUTPUT);
  pinMode(LED_BUILTIN, OUTPUT);
}

void loop() {

  //set LE to low
  digitalWrite(latch, LOW);

  //main code for sending data
  if (count<6) {

    //send data to serial register
    if (attenSet[count]==0){
      digitalWrite(data, LOW);
    } else{
      digitalWrite(data, HIGH);
    }
    
    //clock rise
    digitalWrite(clk, HIGH);

    //blink LED for diagnostics. Fast blinking indicates attenuator is being set.
    digitalWrite(LED_BUILTIN, HIGH);

    //set pulse length
    delay(del);

    //clock fall
    digitalWrite(clk, LOW);
    digitalWrite(LED_BUILTIN, LOW);
    delay(del);
  }

  //toggle LE
  if (count==5) {
    delay(del);
    digitalWrite(latch,HIGH);
    delay(del);
  }
  
  //don't let count get too high
  if (count>5) {
    
    count=10;

    //slow blink LED to signal that loop is complete
    digitalWrite(LED_BUILTIN, HIGH);
    delay(1000);
    digitalWrite(LED_BUILTIN, LOW);
    delay(1000);
  }

  //increment count
  count+=1;
}
