/*=========================================================================================
  ****************************************************************************************
              /=========================================================================/
             /programme:	sub_attack				                                   /
            /filename:		assignment_1_source.cpp                                   /
           /assignment no:	1													     /
		  /version:			0.0001                                                  /
         /hackerjack:		scott ticknor										   /
	    /Student ID:		040431741											  /
       /instructor:			Shawn Unger                                          /
      /course:				cst8132 C++ Programming Language                    /
     /purpose:				to provide a very simple simulation of an epic	   /
    /						sea battle! you are a sub commander whose orders  /
   /							are to eradicate merchant convoy vessels!	 /
  /due by date: 22 Mai 2002												    /
 /=========================================================================/
*****************************************************************************************
=======================================================================================*/




#include <iostream>
#include <conio.h>		//included for getch() only! surely there's a better way...
#include <stdlib.h>
#include <stdio.h>
using namespace std;



//function prototypes...
void Do_Damage(int&, int&, int, int&, int&, int&, int Hit[][100], int Miss[][100], int Summary[][2]);
char Do_Disengage();
char Do_Resupply();
void Do_Mini_Summary(int, int, int Hit[][100], int Miss[][100]);
void Do_Big_Summary(int Summary[][2]);


/*-------------------------------------------------------------------------------------------------
function name:		main
version:			0.1
authour:			Scott Ticknor
studentID			040 431 741
purpose:			program execution :)
					produces random numbers, makes sure counters are
IN from var list:	function prototypes
IN from stdin:		number of torps, choices to disengage/quit
OUT to var list:	??
OUT to stdout:		messages, progress of game, prompts, game data
\-------------------------------------------------------------------------------------------------*/
void main()
{

	//main arrays declared
	int Hit[18][100];
	int Miss[18][100];
	int Summary[18][2];

	//need some main variables
	int torps_i, sunk_i, shots_i, targets_i=6;
	int round_count_i=0, seed=6, flag;  //, hit_var, miss_var;
	char choice_disengage, choice_quit;

	srand(seed);						//does this do anything? only the caribou know...
	targets_i += rand()%7;				//don't know if this is proper usage of rand but it sorta works

	cout<<endl<<"   Welcome aboard U-571, Captain"<<endl;
	cout<<"   Headquarters has issued directives ordering us to "<<endl;
	cout<<"   engage a nearby merchant convoy."<<endl;
	cout<<"   Radar indicates at least 6 target vessels within range, "<<endl;
	cout<<"   and the Helsman has laid in an intercept course..."<<endl;

	cout<<endl<<"   Targets are in sight and weapons are ready. Prepare to Engage!!"<<endl<<endl;

	for (int g=0; g<18; g++)
	{
		Summary[g][0] = 0;
		Summary[g][1] = 0;				//ugly, quick n dirty init here.
	}


	do									//start round
	{

		for (int f=0; f<18; f++)			//preventing garbage values
		{
			for (int g=0; g<100; g++)
			{
				Hit[f][g]  = 0;
				Miss[f][g] = 0;
			}
		}

		torps_i=17;						//re-init these after each re-supply
		sunk_i=0;

			do
			{
			    cout<<"Torpedos remaing: "<<(torps_i+1)<<"  Targets: "<<targets_i<<endl;	//status msg
				cout<<endl<<"Fire how many torpedos? (1,2,3)...>";
				cin>>shots_i;
				while (shots_i!=1 && shots_i!=2 && shots_i!=3)	//no more than 3
				{
					cout<<"Please select 1, 2, or 3 torpedos...>";
					cin>>shots_i;
				}
				while (shots_i > torps_i+1 && torps_i >=0)
				{
					cout<<endl<<"The Ordinance Officer reports we don't have that many torpedoes left!"<<endl;
					cout<<"Please select again"<<endl;
					cin>>shots_i;
				}

				Do_Damage(torps_i, sunk_i, shots_i, flag, targets_i, round_count_i, Hit, Miss, Summary);
				choice_disengage = Do_Disengage();

			}while(choice_disengage=='N' || choice_disengage=='n' && torps_i>=0 && targets_i>0);



   if(targets_i==0)
   {
	  cout<<endl<<"****All targets destroyed! We are victorious!!****"<<endl;
	  Do_Mini_Summary(round_count_i, sunk_i, Hit, Miss); //quick summary after victory
	  cout<<endl;
	  getch();
	  Do_Big_Summary(Summary);
   }

   else if(torps_i<0)
   {
      cout<<endl<<">>> Captain! We are out of torpedos!! <<<"<<endl;
	  cout<<endl<<"Disengaging..."<<endl;
      Do_Mini_Summary(round_count_i, sunk_i, Hit, Miss); //quick summary while empty
   }

   else
   {
      Do_Mini_Summary(round_count_i, sunk_i, Hit, Miss); //quick summary for shits n giggles..oops i mean to meet the requirements
   }




	round_count_i++;				 //increment counter for use in the arrays
	choice_quit = Do_Resupply();

	}while(choice_quit!='q' && choice_quit!='Q'); //if user quits
												  // then print msgs and
	Do_Big_Summary(Summary);					 //  go surfing!
}


//fn Do_Damage
/*-------------------------------------------------------------------------------------------------
function name:		Do_Damage
version:			0.1
authour:			Scott Ticknor
studentID			040 431 741
purpose:			yeeash, this baby does it all!
					loops and calculates Hits and Misses for number of shots_i, stores Hits
					and Misses in respective arrays as well as in array Summary, decrements target
					and torp counters, increments sunk counter
IN from var list:	by reference: torps_i, sunk_i, flag, targets_i, round_count_i,
					arrays Hit, Miss, and Summary
					by value: shots_i
IN from stdin:		none
OUT to var list:	torps_i, sunk_i, flag, targets_i, round_count_i,
					arrays Hit, Miss, and Summary
OUT to stdout:		message IF there is a hit
\-------------------------------------------------------------------------------------------------*/
void Do_Damage(int& torps_i, int& sunk_i, int shots_i, int& flag, int& targets_i,
			   int& round_count_i, int Hit[][100], int Miss[][100], int Summary[][2])
{
   if (torps_i >= 0)
   {
	  flag=0;


      for (int d=0 ; d<shots_i; d++)
	  {

		int torp_enroute=rand()%100;	//this random stuff works like crap
		if (torp_enroute<=33)
		{
			Hit[torps_i][round_count_i] = 1;
			Miss[torps_i][round_count_i] = 0;
   			cout<<"Hit! ";
			flag++;
			Summary[torps_i][0] += 1;		// y-index 0 holds hit values
			torps_i--;
		}
		else
		{
			Hit[torps_i][round_count_i] = 0;
			Miss[torps_i][round_count_i] = 1;
			cout<<"Miss! ";
			Summary[torps_i][1] += 1;		// y-index 1 holds the misses
			torps_i--;
		}
	  }
   }

   if (flag>0)
   {
      cout<<endl<<endl<<"The Tactical Officer reports that one target has been destroyed!"<<endl;
	  cout<<"Well done skipper!"<<endl;
	  sunk_i++;
	  targets_i--;
   }
}



//fn Do_Disengage
/*-------------------------------------------------------------------------------------------------
function name:		Do_Disengage
version:			0.1
authour:			Scott Ticknor
studentID			040 431 741
purpose:			provides user input choice to disengage from attack loop
IN from var list:	none
IN from stdin:		choice (y or n) for variable local_choice
OUT to var list:	returns 'local_choice' to main variable 'choice_disengage'
OUT to stdout:		messages, prompt to choose
\-------------------------------------------------------------------------------------------------*/
char Do_Disengage()
{
	char local_choice;

	cout<<endl<<"Captain, are you ready to disengage (y,n)?...>";
	cin>>local_choice;
	while(local_choice !='y' && local_choice!='Y' && local_choice!='n' && local_choice!='N')
	{
		cout<<"Captain, please tell the Helmsman what to do"<<endl;
		cout<<"disengage (y,n)?...>";
		cin>>local_choice;
	}
	return local_choice;
}


//fn Do_Resupply
/*-------------------------------------------------------------------------------------------------
function name:		Do_Resupply
version:			0.1
authour:			Scott Ticknor
studentID			040 431 741
purpose:			provides user input for making a choice, either quit or restart program execution
IN from var list:	none
IN from stdin:		choice, variable 'local_choice3'
OUT to var list:	function returns 'local_choice3' to main variable 'choice_resupply'
OUT to stdout:		messages, prompt to choose
\-------------------------------------------------------------------------------------------------*/
char Do_Resupply()
{
	char local_choice3;
	cout<<endl<<"Captain, shall we Retreat & Resupply, or Quit and go Surfing (r,q)?...>"<<endl;
	cin>>local_choice3;
	while(local_choice3!='R' && local_choice3!='r' && local_choice3!='q' && local_choice3!='Q')
	{
		cout<<"Captain, we need good judgement at this crucial point in the battle!"<<endl;
		cout<<"What do we do now?"<<endl;
		cout<<"retreat & resupply, or quit (r,q)?...>";
		cin>>local_choice3;
	}
	return local_choice3;
}


//fn Do_Mini_Summary
/*-------------------------------------------------------------------------------------------------
function name:		Do_Mini_Summary
version:			0.1
authour:			Scott Ticknor
studentID			040 431 741
purpose:			reads from Hits and Miss arrays, prints a recap of attack run
IN from var list:	round_count_i, torps_i, sunk_i, int array Hit, int array Miss
IN from stdin:		none
OUT to var list:	none
OUT to stdout:		messages, Hits, Misses, sunk ships, and remaining torps
\-------------------------------------------------------------------------------------------------*/
void Do_Mini_Summary(int round_count_i, int sunk_i, int Hit[][100], int Miss[][100])
{
	int hit_var=0;
	int miss_var=0;
	for (int b=17; b>=0; b--)			//check this shit out..it fixed the garbage values!!!!!
	{
		hit_var  += Hit[b][round_count_i];
		miss_var += Miss[b][round_count_i];
	}
	cout<<endl<<"<---Quick Battle Summary--->"<<endl;
	cout<<"Hits:   "<<hit_var<<endl<<"Misses: "<<miss_var<<endl;
	cout<<"Ships downed: "<<sunk_i<<endl;
}


//fn Do_Big_Summary
/*--------------------------------------------------------------------------------------------------
function name:		Do_Big_Summary
version:			0.1
authour:			Scott Ticknor
student ID:			040 431 741
purpose:			reads values from summary array, prints final results from program procedure
IN from var list:	int array 'Summary'
IN from stdin:		none
OUT to var list:	none
OUT to stdout:		messages, totals from array data
\--------------------------------------------------------------------------------------------------*/
void Do_Big_Summary(int Summary[][2])
{
	int total_hit=0;
	int total_miss=0;

	cout<<"           ----======[ Torpedo Analysis ]======----"<<endl<<endl;

	for (int c=0; c<18; c++) //~:>
	{
		cout<<"Torpedo "<<c<<" Hits: "<<Summary[c][0]<<" Misses: "<<Summary[c][1]<<endl;
	}
	cout<<endl<<"$$ Make Love Not War, man $$"<<endl<<endl;
}

//////////////////////////////////[ E-N-D-=^=-O-F-=^=-F-I-L-E ]////////////////////////////////////////
