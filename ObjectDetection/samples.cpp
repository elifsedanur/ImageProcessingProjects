#include <stdio.h>
#include <stdlib.h>
#include "sqlite3.h" 
#include <iostream>

using namespace std;
sqlite3* db;
sqlite3_stmt* stmt;
int result;

void connection();

int main() {
	connection();
}
void connection() {
	if (sqlite3_open("test.db", &db) == SQLITE_OK) {
		result = sqlite3_prepare_v2(db, " CREATE TABLE IF NOT EXISTS moments (objectName varchar(80),F1	float, F2 float,F3	float, F4	float,F5 float,F6	float, F7 float);", -1, & stmt, NULL);
			sqlite3_step(stmt);
		sqlite3_finalize(stmt);

		if (result != SQLITE_OK) {
			cout << "Error: " << sqlite3_errmsg(db) << "\n";
		}
		else {
			cout << "Table created successfully";
		}
	}
}