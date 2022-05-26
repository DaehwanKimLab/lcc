To compile .proto files run

When inside the /protos directory:

python -m grpc_tools.protoc -I. --python_out=. --grpc_python_out=. ./lccsimulation.proto