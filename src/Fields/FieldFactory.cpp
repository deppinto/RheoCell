#include "FieldFactory.h"
#include "MultiPhaseField.h"

FieldPtr FieldFactory::make_field() {
	// the default box is the cubic one
	char field_type[512] = "multiphasefield";

	if(!strncmp(field_type, "multiphasefield", 512)) throw("Unsupported field", field_type);
	else return std::make_shared<MultiPhaseField>();
}
