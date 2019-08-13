// declare variables
const int LE1 = 2;
const int LE2 = 8;
const int clk1 = 3;
const int clk2 = 9;
const int data1 = 4;
const int data2 = 10;
int count = 0;
int attenSet[6] = {1,0,0,0,0,0};
int del = 50;

void setup() {

  //setup pins
  pinMode(LE1,OUTPUT);
  pinMode(LE2,OUTPUT);
  pinMode(clk1,OUTPUT);
  pinMode(clk2,OUTPUT);
  pinMode(data1,OUTPUT);
  pinMode(data2,OUTPUT);
  pinMode(LED_BUILTIN, OUTPUT);
}

void loop() {

  //set LE to low
  digitalWrite(LE1, LOW);
  digitalWrite(LE2, LOW);

  //main code for sending data
  if (count<6) {

    //send data to serial register
    if (attenSet[count]==0){
      digitalWrite(data1, LOW);
      digitalWrite(data2, LOW);
    } else{
      digitalWrite(data1, HIGH);
      digitalWrite(data2, HIGH);
    }
    
    //clock rise
    digitalWrite(clk1, HIGH);
    digitalWrite(clk2, HIGH);

    //blink LED for diagnostics
    digitalWrite(LED_BUILTIN, HIGH);

    //set pulse length
    delay(del);

    //clock fall
    digitalWrite(clk1, LOW);
    digitalWrite(clk2, LOW);
    digitalWrite(LED_BUILTIN, LOW);
    delay(del);
  }

  //toggle LE
  if (count==5) {
    delay(del);
    digitalWrite(LE1,HIGH);
    digitalWrite(LE2,HIGH);
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
